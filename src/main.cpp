#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/direction_fields.h"
#include "geometrycentral/surface/stripe_patterns.h"
#include "geometrycentral/surface/vector_heat_method.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include <TinyAD/Support/GeometryCentral.hh>
#include "TinyAD/ScalarFunction.hh"
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>

#include "polyscope/pick.h"

#include "geometry.hh"
#include "bspline.hpp"

#include "args/args.hxx"
#include "imgui.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
//std::unique_ptr<ManifoldSurfaceMesh> mesh;
//std::unique_ptr<VertexPositionGeometry> geometry;
std::unique_ptr<ManifoldSurfaceMesh> mesh_bez_surf;
std::unique_ptr<VertexPositionGeometry> geometry_bez_surf;
std::unique_ptr<ManifoldSurfaceMesh> mesh_bez_cnet;
std::unique_ptr<VertexPositionGeometry> geometry_bez_cnet;

// Polyscope visualization handle, to quickly add data to the surface
//polyscope::SurfaceMesh *psMesh;
polyscope::SurfaceMesh* psMesh_bez_surf;
polyscope::SurfaceMesh* psMesh_bez_cnet;
polyscope::CurveNetwork* psCurvNet_bez_cnet;
polyscope::CurveNetwork* psCurvNet_bez_smoothed_cnet;

VertexData<double> us;
VertexData<double> vs;

VertexData<Vector3> dus;
VertexData<Vector3> dvs;
VertexData<Vector3> duus;
VertexData<Vector3> duvs;
VertexData<Vector3> dvvs;

VertexData<Vector3> duus_smoothed;
VertexData<Vector3> duvs_smoothed;
VertexData<Vector3> dvvs_smoothed;

VertexData<double> mean;
VertexData<double> gauss;
VertexData<double> DJ;
VertexData<Vector3> dmin;
VertexData<Vector3> dmax;
VertexData<Vector2> Q;

VertexData<Eigen::Matrix3d> W;
VertexData<double> W_eval1;
VertexData<double> W_eval2;
VertexData<double> W_eval3;
VertexData<Vector3> W_evec1;
VertexData<Vector3> W_evec2;
VertexData<Vector3> W_evec3;


VertexData<Eigen::Matrix3d> W_smoothed;
VertexData<double> W_eval1_smoothed;
VertexData<double> W_eval2_smoothed;
VertexData<double> W_eval3_smoothed;
VertexData<Vector3> W_evec1_smoothed;
VertexData<Vector3> W_evec2_smoothed;
VertexData<Vector3> W_evec3_smoothed;

VertexData<double> mean_smoothed;

VertexData<Vector3> Ns;

EdgeData<double> cot_re;
EdgeData<double> cot_im;


Geometry::BSSurface surf;

// Some algorithm parameters
float param1 = 42.0;
int num_smoothing_iters = 0;
float w_ders = 1.0;
float w_oldcps = 1.0;

Vector3 inSystem(const Vector3& u, const Vector3& v, const Vector3& w) {
    // Algebraic version
    double v2 = dot(v, v), w2 = dot(w, w), uv = dot(u , v), uw = dot(u, w), vw = dot(v, w);
    double denom = v2 * w2 - vw * vw;
    return { (w2 * uv - vw * uw)/denom, (v2 * uw - vw * uv)/denom, 0 };
    // Geometric version (should be the same)
    // auto uv = u ^ v, uw = u ^ w, vw = v ^ w;
    // return {
    //   uw.norm() / vw.norm() * (uw * vw > 0 ? 1 : -1),
    //   uv.norm() / vw.norm() * (uv * vw < 0 ? 1 : -1),
    //   0
    // };
}

SparseMatrix<std::complex<double>> computeVertexConnectionLaplacian(IntrinsicGeometryInterface& geometry, int nSym) {

    SurfaceMesh& mesh = geometry.mesh;

    geometry.requireVertexIndices();
    geometry.requireEdgeCotanWeights();
    geometry.requireTransportVectorsAlongHalfedge();

    cot_re = EdgeData<double>(*mesh_bez_surf);
    cot_im = EdgeData<double>(*mesh_bez_surf);

    std::vector<Eigen::Triplet<std::complex<double>>> triplets;
    for (Halfedge he : mesh.halfedges()) {

        size_t iTail = geometry.vertexIndices[he.vertex()];
        size_t iTip = geometry.vertexIndices[he.next().vertex()];

        // Levi-Civita connection between vertices
        Vector2 rot = geometry.transportVectorsAlongHalfedge[he.twin()].pow(nSym);
        double weight = geometry.edgeCotanWeights[he.edge()];

        cot_re[he.edge()] = (weight * rot)[0];
        cot_im[he.edge()] = (weight * rot)[1];
        triplets.emplace_back(iTail, iTail, weight);
        triplets.emplace_back(iTail, iTip, -weight * rot);
    }
    psMesh_bez_surf->addEdgeScalarQuantity("cot (Re)", cot_re);
    psMesh_bez_surf->addEdgeScalarQuantity("cot (Im)", cot_im);
    // assemble matrix from triplets
    Eigen::SparseMatrix<std::complex<double>> vertexConnectionLaplacian(mesh.nVertices(), mesh.nVertices());
    vertexConnectionLaplacian.setFromTriplets(triplets.begin(), triplets.end());

    // Shift to avoid singularity
    Eigen::SparseMatrix<std::complex<double>> eye(mesh.nVertices(), mesh.nVertices());
    eye.setIdentity();
    vertexConnectionLaplacian += 1e-9 * eye;

    return vertexConnectionLaplacian;
}

void fitBeziertoDerivatives(
    const std::vector<Eigen::Vector3d> &der_values, 
    const std::vector<Eigen::Vector2d> &der_uvpars,
    size_t                             deg_u,
    size_t                             deg_v,
    const std::vector<Eigen::Vector3d> &old_cps,
    std::vector<Eigen::Vector3d>       &new_cps,
    size_t                             dd          = 2,
    double                             w_ders_      = 1.0,
    double                             w_oldcps_    = 1.0
)
{
    size_t ncp = old_cps.size();
    size_t nd = der_values.size()/3;
    SparseMatrix<double> MM(3*nd + ncp, ncp);
    //Eigen::VectorXi MM_nnz = Eigen::VectorXi::Constant(3 * nd + ncp, ncp);
    //for (size_t i = 0; i < ncp; ++i) {
    //    MM_nnz(3 * nd + i) = 1;
    //}
    //MM.reserve(MM_nnz);

    SparseMatrix<double> WW(3*nd + ncp, 3*nd + ncp);
    WW.reserve(Eigen::VectorXd::Constant(3 * nd + ncp, 1));
    DenseMatrix<double> rhs(3*nd + ncp, 3);
    BSpline::PointVector zero_cps = BSpline::zeroPoints(BSpline::fromEigen(old_cps));
    for (size_t i = 0; i < nd; ++i) {
        double u = der_uvpars[3 * i][0];
        double v = der_uvpars[3 * i][1];
        for (size_t ci = 0; ci < ncp; ++ci) {
            auto cps4der = zero_cps;
            cps4der[ci] = { 1.0, 0.0, 0.0 };
            BSpline::BSSurface surf_temp(deg_u, deg_v, cps4der);
            BSpline::PointMatrix der;
            surf_temp.eval(u, v, dd, der);
            const auto& duu = der[2][0];
            const auto& duv = der[1][1];
            const auto& dvv = der[0][2];
            MM.coeffRef(3 * i, ci)     = duu[0];
            MM.coeffRef(3 * i + 1, ci) = duv[0];
            MM.coeffRef(3 * i + 2, ci) = dvv[0];
        }
        WW.coeffRef(3 * i,     3 * i)     = w_ders_;
        WW.coeffRef(3 * i + 1, 3 * i + 1) = w_ders_;
        WW.coeffRef(3 * i + 2, 3 * i + 2) = w_ders_;
        rhs.row(3 * i)     = der_values[3 * i];
        rhs.row(3 * i + 1) = der_values[3 * i + 1];
        rhs.row(3 * i + 2) = der_values[3 * i + 2];

    }
    for (size_t ci = 0; ci < ncp; ++ci) {
        MM.coeffRef(3 * nd + ci, ci)          = 1.0;
        WW.coeffRef(3 * nd + ci, 3 * nd + ci) = w_oldcps_;
        rhs.row(3 * nd + ci) = old_cps[ci];
    }
    Eigen::SparseMatrix<double> QQ = MM.transpose() * WW * MM;
    Eigen::SparseLU solver(QQ);
    DenseMatrix<double> x = solver.solve(MM.transpose()*WW*rhs);
    std::cerr << "Residual (der): " << (MM.topRows(3*nd) * x.topRows(3*nd) - rhs.topRows(3*nd)).col(2).norm() << std::endl;
    std::cerr << "Residual (CPs): " << (MM.bottomRows(ncp) * x.bottomRows(ncp) - rhs.bottomRows(ncp)).col(2).norm() << std::endl;
    new_cps = old_cps;
    for (size_t ci = 0; ci < ncp; ++ci) {
        new_cps[ci] = { x(ci,0), x(ci,1), x(ci,2) };
    }
}

// Example computation function -- this one computes and registers a scalar
// quantity
void doWork() {
  //polyscope::warning("Computing Gaussian curvature.\nalso, parameter value = " +
  //                   std::to_string(param1));


    //geometry->requireEdgeCotanWeights();
    //auto func = TinyAD::scalar_function<3>(mesh->vertices());

    //// Add an objective term per triangle. Each connecting 3 vertices
    //func.add_elements<2>(mesh->edges(), [&](auto& element)->TINYAD_SCALAR_TYPE(element)
    //{
    //    // Element is evaluated with either double or TinyAD::Double<6>
    //    using T = TINYAD_SCALAR_TYPE(element);

    //    Edge e = element.handle;

    //    //Vertex v = element.handle;
    //    Eigen::Vector3<T> p0 = element.variables(e.firstVertex());
    //    Eigen::Vector3<T> p1 = element.variables(e.secondVertex());

    //    return 0.5 * geometry->edgeCotanWeights[e] * (p1 - p0).squaredNorm();
    //}
    //);

    //Eigen::VectorXd x = func.x_from_data([&](Vertex v)->Eigen::Vector3d {
    //    auto p = geometry->vertexPositions[v];
    //    return { p[0], p[1], p[2] };
    //    });


    //// Projected Newton
    //TinyAD::LinearSolver solver;
    //int max_iters = 10;
    //double convergence_eps = 1e-6;
    //for (int i = 0; i < max_iters; ++i)
    //{
    //    auto [f, g] = func.eval_with_gradient(x);
    //    TINYAD_DEBUG_OUT("Energy before " << i << ": " << f);
    //    double lr = 1.0;
    //    Eigen::VectorXd d = -lr * g;//TinyAD::newton_direction(g, H_proj, solver, 0.01);
    //    //if (TinyAD::newton_decrement(d, g) < convergence_eps)
    //    //    break;
    //    x = TinyAD::line_search(x, d, f, g, func);

    //    //VertexData<Vector3> grad(*mesh);
    //    //func.x_to_data(g, [&](Vertex v, const Eigen::Vector3d& p) {
    //    //    grad[v] = { p[0], p[1], p[2] };
    //    //    }
    //    //);

    //    func.x_to_data(x, [&](Vertex v, const Eigen::Vector3d& p) {
    //        geometry->vertexPositions[v] = { p[0], p[1], p[2] };
    //        }
    //    );
    //    geometry->refreshQuantities();
    //    //psMesh->updateVertexPositions(geometry->vertexPositions);

    //    //psMesh->addVertexVectorQuantity("Gradient", grad);
    //}
    //TINYAD_DEBUG_OUT("Energy after: " << func.eval(x));

    //func.x_to_data(x, [&](Vertex v, const Eigen::Vector3d& p) {
    //    geometry->vertexPositions[v] = { p[0], p[1], p[2] };
    //    }
    //);
    //psMesh->updateVertexPositions(geometry->vertexPositions);




    Geometry::PointVector cps;
    std::vector<Vector3> cps_gc;
    size_t deg_u = 5;
    size_t deg_v = 5;
    for (int i = 0; i <= deg_u; ++i) {
        for (int j = 0; j <= deg_v; ++j) {
            double x = (double(i) - double(deg_u)/2)/(double(deg_u) / 2);
            double y = (double(j) - double(deg_v)/2)/(double(deg_v) / 2);
            Vector3 cp = { x, y, randomReal(-1, 1) /*x * x + y * y + randomReal(-1.0, 0.0)*/ };
            cps.push_back({ cp[0], cp[1], cp[2] });
            cps_gc.push_back(cp);
        }
    }

    surf = Geometry::BSSurface(deg_u, deg_v, cps);

    size_t resolution = 100;
    BSpline::writeToMesh(surf, resolution);

    std::tie(mesh_bez_surf, geometry_bez_surf) = readManifoldSurfaceMesh("bezier.obj");

    psMesh_bez_surf = polyscope::registerSurfaceMesh(
        polyscope::guessNiceNameFromPath("bezier.obj"),
        geometry_bez_surf->inputVertexPositions, mesh_bez_surf->getFaceVertexList());

    psMesh_bez_surf->setAllPermutations(polyscopePermutations(*mesh_bez_surf));


    // Add control net as mesh
    //std::vector<std::vector<size_t>> polygons;
    //for (size_t i = 0; i < deg_u; ++i) {
    //    for (size_t j = 0; j < deg_v; ++j) {
    //        polygons.push_back({
    //        i * (deg_u + 1) + j,
    //        i * (deg_u + 1) + j + 1,
    //        (i + 1) * (deg_u + 1) + j + 1,
    //        (i + 1) * (deg_u + 1) + j });
    //    }
    //}

    //std::tie(mesh_bez_cnet, geometry_bez_cnet) = geometrycentral::surface::makeManifoldSurfaceMeshAndGeometry(polygons, cps_gc);

    //psMesh_bez_cnet = polyscope::registerSurfaceMesh(
    //    "Control Polyhedron",
    //    geometry_bez_cnet->inputVertexPositions, mesh_bez_cnet->getFaceVertexList());


    // Add control net as curve network
    std::vector<std::array<size_t, 2>> cp_edges;
    for (size_t i = 0; i <= deg_u; ++i) {
        for (size_t j = 0; j < deg_v; ++j) {
            cp_edges.push_back({ i * (deg_u + 1) + j , i * (deg_u + 1) + j + 1 });
            if (i < deg_u) {
                cp_edges.push_back({ i * (deg_u + 1) + j ,(i + 1) * (deg_u + 1) + j });
                if (j == deg_v - 1) {
                    cp_edges.push_back({ i * (deg_u + 1) + j + 1 ,(i + 1) * (deg_u + 1) + j + 1 });
                }
            }
        }
    }
    psCurvNet_bez_cnet = polyscope::registerCurveNetwork("Control net", cps_gc, cp_edges);
    psCurvNet_bez_smoothed_cnet = polyscope::registerCurveNetwork("Control net (smoothed)", cps_gc, cp_edges);

    //psMesh_bez_cnet->setAllPermutations(polyscopePermutations(*mesh_bez_cnet));


    //VertexData<double> us(*mesh_bez_surf);
    //VertexData<double> vs(*mesh_bez_surf);

    //VertexData<Vector3> dus(*mesh_bez_surf);
    //VertexData<Vector3> dvs(*mesh_bez_surf);
    //VertexData<Vector3> duus(*mesh_bez_surf);
    //VertexData<Vector3> duvs(*mesh_bez_surf);
    //VertexData<Vector3> dvvs(*mesh_bez_surf);

    //VertexData<double> mean(*mesh_bez_surf);
    //VertexData<double> gauss(*mesh_bez_surf);
    //VertexData<double> DJ(*mesh_bez_surf);
    //VertexData<Vector2> dmin(*mesh_bez_surf);
    //VertexData<Vector2> dmax(*mesh_bez_surf);
    //VertexData<Vector2> Q(*mesh_bez_surf);

    //VertexData<Vector3> Ns(*mesh_bez_surf);

    us = VertexData<double>(*mesh_bez_surf);
    vs = VertexData<double>(*mesh_bez_surf);

    dus  = VertexData<Vector3>(*mesh_bez_surf);
    dvs  = VertexData<Vector3>(*mesh_bez_surf);
    duus = VertexData<Vector3>(*mesh_bez_surf);
    duvs = VertexData<Vector3>(*mesh_bez_surf);
    dvvs = VertexData<Vector3>(*mesh_bez_surf);

    mean  = VertexData<double>(*mesh_bez_surf);
    gauss = VertexData<double>(*mesh_bez_surf);
    DJ    = VertexData<double>(*mesh_bez_surf);
    dmin  = VertexData<Vector3>(*mesh_bez_surf);
    dmax  = VertexData<Vector3>(*mesh_bez_surf);
    Q     = VertexData<Vector2>(*mesh_bez_surf);
    VertexData<Vector3> ev1_Q(*mesh_bez_surf);
    VertexData<Vector3> ev2_Q(*mesh_bez_surf);
    VertexData<double> eval1_Q(*mesh_bez_surf);
    VertexData<double> eval2_Q(*mesh_bez_surf);
    VertexData<double> kmin(*mesh_bez_surf);
    VertexData<double> kmax(*mesh_bez_surf);

    W     = VertexData<Eigen::Matrix3d>(*mesh_bez_surf);
    W_eval1 = VertexData<double>(*mesh_bez_surf);
    W_eval2 = VertexData<double>(*mesh_bez_surf);
    W_eval3 = VertexData<double>(*mesh_bez_surf);
    W_evec1 = VertexData<Vector3>(*mesh_bez_surf);
    W_evec2 = VertexData<Vector3>(*mesh_bez_surf);
    W_evec3 = VertexData<Vector3>(*mesh_bez_surf);

    Ns = VertexData<Vector3>(*mesh_bez_surf);

    VertexData<Vector3> b1(*mesh_bez_surf);
    VertexData<Vector3> b2(*mesh_bez_surf);
    VertexData<Vector2> dmin_W(*mesh_bez_surf);
    VertexData<Vector2> dmax_W(*mesh_bez_surf);

    geometry_bez_surf->requireVertexNormals();
    geometry_bez_surf->requireVertexDualAreas();

    psMesh_bez_surf->vertexNormals.ensureHostBufferAllocated();

    for (size_t i = 0; i < resolution; ++i) {
        double v = (double)i / (double)(resolution - 1);
        for (size_t j = 0; j < resolution; ++j) {
            double u = (double)j / (double)(resolution - 1);
            
            BSpline::VectorMatrix der;
            surf.eval(u, v, 2, der);
            const auto& du = der[1][0];
            const auto& dv = der[0][1];
            const auto& duu = der[2][0];
            const auto& duv = der[1][1];
            const auto& dvv = der[0][2];
            double E = du * du;
            double F = du * dv;
            double G = dv * dv;
            auto n = du ^ dv;
            auto n_ = n.normalized();
            double L = n_ * duu;
            double M = n_ * duv;
            double N = n_ * dvv;

            Eigen::MatrixXd J(3, 2);
            J.col(0) = Eigen::Vector3d( du[0], du[1], du[2] );
            J.col(1) = Eigen::Vector3d( dv[0], dv[1], dv[2] );

            Eigen::Matrix2d I;
            I(0, 0) = E; I(0, 1) = F;
            I(1, 0) = F; I(1, 1) = G;

            Eigen::MatrixXd Jpinv = I.inverse() * J.transpose();

            Eigen::Matrix2d II;
            II(0, 0) = L; II(0, 1) = M;
            II(1, 0) = M; II(1, 1) = N;

            Eigen::Matrix2d S = I.inverse()*II;
            Eigen::Matrix3d W_ = Jpinv.transpose() * II * Jpinv;

            size_t v_idx = i * resolution + j;

            us[v_idx] = u;
            vs[v_idx] = v;

            dus[v_idx] = { du[0], du[1], du[2] };
            dvs[v_idx] = { dv[0], dv[1], dv[2] };

            duus[v_idx] = { duu[0], duu[1], duu[2] };
            duvs[v_idx] = { duv[0], duv[1], duv[2] };
            dvvs[v_idx] = { dvv[0], dvv[1], dvv[2] };

            double H = S.trace() / 2.0;
            double K = S.determinant();

            mean[v_idx] = H;
            gauss[v_idx] = K;

            Ns[v_idx] = { n_[0], n_[1], n_[2] };
            DJ[v_idx] = n.norm();

            auto bs = Ns[v_idx].buildTangentBasis();
            b1[v_idx] = bs[0];
            b2[v_idx] = bs[1];

            
            Eigen::EigenSolver<Eigen::Matrix2d> es_S(S);
            Eigen::Vector2d e1_uv = es_S.eigenvectors().real().col(0);
            Eigen::Vector2d e2_uv = es_S.eigenvectors().real().col(1);
            auto e1_xyz = (du * e1_uv[0] + dv * e1_uv[1]).normalize();
            auto e2_xyz = (du * e2_uv[0] + dv * e2_uv[1]).normalize();
            dmin[v_idx] = { e1_xyz[0], e1_xyz[1], e1_xyz[2] };
            dmax[v_idx] = { e2_xyz[0], e2_xyz[1], e2_xyz[2] };
            //dmin[v_idx] = { e1_uv[0], e1_uv[1] };
            //dmax[v_idx] = { e2_uv[0], e2_uv[1] };
            kmin[v_idx] = es_S.eigenvalues().real()[0];
            kmax[v_idx] = es_S.eigenvalues().real()[1];
            Eigen::Matrix2d Hopf = S - Eigen::Matrix2d::Identity() * H;
            Eigen::EigenSolver<Eigen::Matrix2d> es_Q(Hopf);
            //kmin[v_idx] = std::abs(es_Q.eigenvalues().real()[0]);
            //kmax[v_idx] = std::abs(es_Q.eigenvalues()[1].real());
            eval1_Q[v_idx] = es_Q.eigenvalues().real()[0];
            eval2_Q[v_idx] = es_Q.eigenvalues().real()[1];
            size_t Q1_idx = eval1_Q[v_idx] < eval2_Q[v_idx] ? 0 : 1;
            size_t Q2_idx = eval1_Q[v_idx] < eval2_Q[v_idx] ? 1 : 0;
            eval1_Q[v_idx] = es_Q.eigenvalues().real()[Q1_idx];
            eval2_Q[v_idx] = es_Q.eigenvalues().real()[Q2_idx];
            Eigen::Vector3d e1_Q_xyz = (J*(es_Q.eigenvectors().col(Q1_idx).real())).normalized();
            Eigen::Vector3d e2_Q_xyz = (J*(es_Q.eigenvectors().col(Q2_idx).real())).normalized();
            ev1_Q[v_idx] = { e1_Q_xyz(0), e1_Q_xyz(1), e1_Q_xyz(2) };
            ev2_Q[v_idx] = { e2_Q_xyz(0), e2_Q_xyz(1), e2_Q_xyz(2) };

            Q[v_idx] = { dot(ev1_Q[v_idx], bs[0]), dot(ev1_Q[v_idx], bs[1]) };
            Q[v_idx] *= std::abs(es_Q.eigenvalues().real()[Q1_idx]);

            // Weingarten
            W[v_idx] = W_;
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es_W(W_);

            auto evals = es_W.eigenvalues();
            double eps = 1E-12;
            size_t idx_zero = std::abs(evals[0]) < eps ? 0 : (std::abs(evals[1]) < eps ? 1 : 2);
            size_t idx_k1 = 1;
            size_t idx_k2 = 2;
            if (idx_zero == 1) {
                idx_k1 = 0;
                idx_k2 = 2;
            }
            else if (idx_zero == 2) {
                idx_k1 = 0;
                idx_k2 = 1;
            }


            W_eval1[v_idx] = evals[idx_k1];
            W_eval2[v_idx] = evals[idx_k2];
            W_eval3[v_idx] = evals[idx_zero];
            Eigen::Vector3d e1 = es_W.eigenvectors().col(idx_k1);
            Eigen::Vector3d e2 = es_W.eigenvectors().col(idx_k2);
            Eigen::Vector3d e3 = es_W.eigenvectors().col(idx_zero);
            W_evec1[v_idx] = { e1[0], e1[1], e1[2] };
            W_evec2[v_idx] = { e2[0], e2[1], e2[2] };
            W_evec3[v_idx] = { e3[0], e3[1], e3[2] };
            W_evec3[v_idx] *= geometrycentral::dot(W_evec3[v_idx], Ns[v_idx]);

            dmin_W[v_idx] = { dot(bs[0],  W_evec1[v_idx]), dot(bs[1],  W_evec1[v_idx]) };
            dmax_W[v_idx] = { dot(bs[0],  W_evec2[v_idx]), dot(bs[1],  W_evec2[v_idx]) };

            geometry_bez_surf->vertexNormals[v_idx] = Ns[v_idx];

        }
    }

    psMesh_bez_surf->addVertexScalarQuantity("u", us, polyscope::DataType::MAGNITUDE);
    psMesh_bez_surf->addVertexScalarQuantity("v", vs, polyscope::DataType::MAGNITUDE);
    psMesh_bez_surf->addVertexVectorQuantity("Du", dus);
    psMesh_bez_surf->addVertexVectorQuantity("Dv", dvs);
    //psMesh_bez_surf->addVertexVectorQuantity("Duu", duus);
    //psMesh_bez_surf->addVertexVectorQuantity("Duv", duvs);
    //psMesh_bez_surf->addVertexVectorQuantity("Dvv", dvvs);
    psMesh_bez_surf->addVertexVectorQuantity("Vertex normals (mesh)", geometry_bez_surf->vertexNormals);
    psMesh_bez_surf->addVertexVectorQuantity("Vertex normals (exact)", Ns);
    //psMesh_bez_surf->addVertexTangentVectorQuantity("dmin", dmin, dus, dvs);
    //psMesh_bez_surf->addVertexTangentVectorQuantity("dmax", dmax, dus, dvs);
    psMesh_bez_surf->addVertexVectorQuantity("dmin", dmin);
    psMesh_bez_surf->addVertexVectorQuantity("dmax", dmax);
    psMesh_bez_surf->addVertexScalarQuantity("kmin",
        kmin,
        polyscope::DataType::SYMMETRIC);
    psMesh_bez_surf->addVertexScalarQuantity("kmax",
        kmax,
        polyscope::DataType::SYMMETRIC);
    psMesh_bez_surf->addVertexTangentVectorQuantity("Q", Q, b1, b2, 2);
    psMesh_bez_surf->addVertexVectorQuantity("Q eigenvectors", ev1_Q);
    //psMesh_bez_surf->addVertexScalarQuantity("Q eigenvalues 1",
    //    eval1_Q,
    //    polyscope::DataType::MAGNITUDE);
    //psMesh_bez_surf->addVertexScalarQuantity("Q eigenvalues 2",
    //    eval2_Q,
    //    polyscope::DataType::MAGNITUDE);
    psMesh_bez_surf->addVertexScalarQuantity("Q magnitude",
        eval2_Q,
        polyscope::DataType::MAGNITUDE);

    psMesh_bez_surf->addVertexVectorQuantity("W eigenvectors 2", W_evec2);
    psMesh_bez_surf->addVertexVectorQuantity("W eigenvectors 1", W_evec1);
    psMesh_bez_surf->addVertexVectorQuantity("W eigenvectors 2", W_evec2);
    psMesh_bez_surf->addVertexTangentVectorQuantity("W eigenvectors 1 (local)", dmin_W, b1, b2, 2);
    psMesh_bez_surf->addVertexTangentVectorQuantity("W eigenvectors 2 (local)", dmax_W, b1, b2, 2);
    psMesh_bez_surf->addVertexVectorQuantity("W eigenvectors 3", W_evec3);
    psMesh_bez_surf->addVertexScalarQuantity("W eigenvalues 1",
        W_eval1,
        polyscope::DataType::SYMMETRIC);
    psMesh_bez_surf->addVertexScalarQuantity("W eigenvalues 2",
        W_eval2,
        polyscope::DataType::SYMMETRIC);
    psMesh_bez_surf->addVertexScalarQuantity("W eigenvalues 3",
        W_eval3,
        polyscope::DataType::SYMMETRIC);
    //static_cast<polyscope::SurfaceVertexTangentVectorQuantity*>(psMesh->getQuantity("Q"))->setVectorRadius(
    //    static_cast<polyscope::SurfaceVertexTangentVectorQuantity*>(psMesh->getQuantity("Q"))->getVectorRadius(), false
    //);

    //geometry->requireVertexGaussianCurvatures();
    //geometry->requireVertexMeanCurvatures();
    //geometry->requireEdgeLengths();
    psMesh_bez_surf->addVertexScalarQuantity("Gaussian curvature",
        gauss,
        polyscope::DataType::SYMMETRIC);
    psMesh_bez_surf->addVertexScalarQuantity("Mean curvature",
        mean,
        polyscope::DataType::SYMMETRIC);
    //psMesh_bez_surf->addVertexScalarQuantity("DJ",
    //    DJ,
    //    polyscope::DataType::MAGNITUDE);
    //psMesh->addEdgeScalarQuantity("edgelength",
    //    geometry->edgeLengths,
    //    polyscope::DataType::MAGNITUDE);

    //auto picked = polyscope::pick::getSelection();
    //if (picked.first != nullptr) {
    //}
    psMesh_bez_surf->setShadeStyle(polyscope::MeshShadeStyle::Smooth);
    psMesh_bez_surf->refresh();
    polyscope::view::resetCameraToHomeView();
}

void refreshNormals() {

    for (auto v : mesh_bez_surf->vertices()) {
        size_t v_idx = v.getIndex();
        auto Ns = static_cast<polyscope::SurfaceVertexVectorQuantity*>(psMesh_bez_surf->getQuantity("Vertex normals (exact)"));
        psMesh_bez_surf->vertexNormals.data[v_idx] = Ns->vectors.data[v_idx] ;
        geometry_bez_surf->vertexNormals[v_idx] = { Ns->vectors.data[v_idx][0], Ns->vectors.data[v_idx][1], Ns->vectors.data[v_idx][2] };

    }
    psMesh_bez_surf->vertexNormals.markHostBufferUpdated();
    //psMesh_bez_surf->refresh();
}

void smoothDerivatives()
{
    //auto duus = VertexData<Vector3>(mesh_bez_surf, psMesh_bez_surf->getQuantity("Duu"));
    duus_smoothed = VertexData<Vector3>(*mesh_bez_surf);
    duvs_smoothed = VertexData<Vector3>(*mesh_bez_surf);
    dvvs_smoothed = VertexData<Vector3>(*mesh_bez_surf);
    VectorHeatMethodSolver vhs(*geometry_bez_surf);
    geometry_bez_surf->requireCotanLaplacian();
    geometry_bez_surf->requireVertexLumpedMassMatrix();
    const auto& LL = geometry_bez_surf->cotanLaplacian;
    const auto& MM = geometry_bez_surf->vertexLumpedMassMatrix;
    double dt = 1.0;
    SparseMatrix<double> heat = MM + dt * LL;
    PositiveDefiniteSolver solv(heat);
    for (size_t cc = 0; cc < 3; ++cc) {
        VertexData<double> rhs_uu(*mesh_bez_surf);
        VertexData<double> rhs_uv(*mesh_bez_surf);
        VertexData<double> rhs_vv(*mesh_bez_surf);
        for (auto v : mesh_bez_surf->vertices()) {
            rhs_uu[v] = duus[v][cc];
            rhs_uv[v] = duvs[v][cc];
            rhs_vv[v] = dvvs[v][cc];
        }
        //auto MM = geometry_bez_surf->vertexLumpedMassMatrix;
        //
        //VertexData<double> x_pw(*mesh_bez_surf, solv.solve(x.toVector()));
        //VertexData<double> x(*mesh_bez_surf, solv.solve(MM*(rhs.toVector())));
        auto x_uu = rhs_uu;
        auto x_uv = rhs_uv;
        auto x_vv = rhs_vv;
        int num_iter = num_smoothing_iters;
        while (num_iter-- > 0) {
            x_uu = vhs.scalarDiffuse(VertexData<double>(*mesh_bez_surf, MM * x_uu.toVector()));
            x_uv = vhs.scalarDiffuse(VertexData<double>(*mesh_bez_surf, MM * x_uv.toVector()));
            x_vv = vhs.scalarDiffuse(VertexData<double>(*mesh_bez_surf, MM * x_vv.toVector()));
        }

        //VertexData<double> x = vhs.scalarDiffuse());
        for (auto v : mesh_bez_surf->vertices()) {
            duus_smoothed[v][cc] = x_uu[v];
            duvs_smoothed[v][cc] = x_uv[v];
            dvvs_smoothed[v][cc] = x_vv[v];
        }
    }
    psMesh_bez_surf->addVertexVectorQuantity("Duu (smoothed)", duus_smoothed);
    psMesh_bez_surf->addVertexVectorQuantity("Duv (smoothed)", duvs_smoothed);
    psMesh_bez_surf->addVertexVectorQuantity("Dvv (smoothed)", dvvs_smoothed);
}

void smoothWeingartens()
{
    W_smoothed = VertexData<Eigen::Matrix3d>(*mesh_bez_surf);
    W_eval1_smoothed = VertexData<double>(*mesh_bez_surf);
    W_eval2_smoothed = VertexData<double>(*mesh_bez_surf);
    W_eval3_smoothed = VertexData<double>(*mesh_bez_surf);
    W_evec1_smoothed = VertexData<Vector3>(*mesh_bez_surf);
    W_evec2_smoothed = VertexData<Vector3>(*mesh_bez_surf);
    W_evec3_smoothed = VertexData<Vector3>(*mesh_bez_surf);

    VertexData<Vector3> basis1(*mesh_bez_surf);
    VertexData<Vector3> basis2(*mesh_bez_surf);
    VertexData<Vector2> dmin_W(*mesh_bez_surf);
    VertexData<Vector2> dmax_W(*mesh_bez_surf);

    geometry_bez_surf->requireVertexLumpedMassMatrix();
    const auto& MM = geometry_bez_surf->vertexLumpedMassMatrix;

    VectorHeatMethodSolver vhs(*geometry_bez_surf);

    for (size_t row = 0; row < 3; ++row) {
        for (size_t col = 0; col < 3; ++col) {
            VertexData<double> rhs(*mesh_bez_surf);
            for (auto v : mesh_bez_surf->vertices()) {
                rhs[v.getIndex()] = W[v.getIndex()](row, col);
            }
            auto x = rhs;
            int num_iter = num_smoothing_iters;
            while (num_iter-- > 0) {
                x = vhs.scalarDiffuse(VertexData<double>(*mesh_bez_surf, MM * x.toVector()));
            }
            for (auto v : mesh_bez_surf->vertices()) {
                W_smoothed[v.getIndex()](row, col) = x[v.getIndex()];
            }
        }
    }

    for (auto v : mesh_bez_surf->vertices()) {

        const auto &W_ = W_smoothed[v];
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es_W(W_);

        auto evals = es_W.eigenvalues();
        double eps = 1E-12;
        size_t idx_zero = std::abs(evals[0]) < eps ? 0 : (std::abs(evals[1]) < eps ? 1 : 2);
        size_t idx_k1 = 1;
        size_t idx_k2 = 2;
        if (idx_zero == 1) {
            idx_k1 = 0;
            idx_k2 = 2;
        }
        else if (idx_zero == 2) {
            idx_k1 = 0;
            idx_k2 = 1;
        }


        W_eval1_smoothed[v] = evals[idx_k1];
        W_eval2_smoothed[v] = evals[idx_k2];
        W_eval3_smoothed[v] = evals[idx_zero];
        Eigen::Vector3d e1 = es_W.eigenvectors().col(idx_k1);
        Eigen::Vector3d e2 = es_W.eigenvectors().col(idx_k2);
        Eigen::Vector3d e3 = es_W.eigenvectors().col(idx_zero);
        W_evec1_smoothed[v] = { e1[0], e1[1], e1[2] };
        W_evec2_smoothed[v] = { e2[0], e2[1], e2[2] };
        W_evec3_smoothed[v] = { e3[0], e3[1], e3[2] };
        //W_evec3_smoothed[v] *= geometrycentral::dot(W_evec3_smoothed[v], Ns[v]);
        auto bs = W_evec3[v].buildTangentBasis();
        basis1[v] = bs[0];
        basis2[v] = bs[1];

        dmin_W[v] = { dot(bs[0],  W_evec1_smoothed[v]), dot(bs[1],  W_evec1_smoothed[v]) };
        dmax_W[v] = { dot(bs[0],  W_evec2_smoothed[v]), dot(bs[1],  W_evec2_smoothed[v]) };
    }

    psMesh_bez_surf->addVertexVectorQuantity("W eigenvectors 1 (smoothed)", W_evec1_smoothed);
    psMesh_bez_surf->addVertexVectorQuantity("W eigenvectors 2 (smoothed)", W_evec2_smoothed);
    psMesh_bez_surf->addVertexTangentVectorQuantity("W eigenvectors 1 (local) (smoothed)", dmin_W, basis1, basis2, 2);
    psMesh_bez_surf->addVertexTangentVectorQuantity("W eigenvectors 2 (local) (smoothed)", dmax_W, basis1, basis2, 2);
    psMesh_bez_surf->addVertexVectorQuantity("W eigenvectors 3 (smoothed)", W_evec3_smoothed);
    psMesh_bez_surf->addVertexScalarQuantity("W eigenvalues 1 (smoothed)",
        W_eval1_smoothed,
        polyscope::DataType::SYMMETRIC);
    psMesh_bez_surf->addVertexScalarQuantity("W eigenvalues 2 (smoothed)",
        W_eval2_smoothed,
        polyscope::DataType::SYMMETRIC);
    psMesh_bez_surf->addVertexScalarQuantity("W eigenvalues 3 (smoothed)",
        W_eval3_smoothed,
        polyscope::DataType::SYMMETRIC);
}

void smoothQs()
{
    duus_smoothed = VertexData<Vector3>(*mesh_bez_surf);
    duvs_smoothed = VertexData<Vector3>(*mesh_bez_surf);
    dvvs_smoothed = VertexData<Vector3>(*mesh_bez_surf);

    geometry_bez_surf->requireEdgeLengths();
    geometry_bez_surf->requireVertexLumpedMassMatrix();
    geometry_bez_surf->requireVertexGalerkinMassMatrix();
    SparseMatrix<Eigen::dcomplex> MM = geometry_bez_surf->vertexLumpedMassMatrix.cast<std::complex<double>>();
    //SparseMatrix<Eigen::dcomplex> MM = geometry_bez_surf->vertexGalerkinMassMatrix.cast<std::complex<double>>();
    auto InvM = MM;
    auto LL = computeVertexConnectionLaplacian(*geometry_bez_surf, 2);
    // Compute mean edge length and set shortTime
    double tCoef = 1.0;
    double meanEdgeLength = 0.;
    for (Edge e : mesh_bez_surf->edges()) {
        meanEdgeLength += geometry_bez_surf->edgeLengths[e];
    }
    meanEdgeLength /= mesh_bez_surf->nEdges();
    double shortTime = tCoef* meanEdgeLength* meanEdgeLength;

    for (auto v : mesh_bez_surf->vertices()) {
        //MM.coeffRef(v.getIndex(), v.getIndex()) = v.degree();
        auto Mii = MM.coeffRef(v.getIndex(), v.getIndex());
        InvM.coeffRef(v.getIndex(), v.getIndex()) = 1.0 / Mii;
    }
    SparseMatrix<Eigen::dcomplex> id = MM;
    id.setIdentity();

    // Build the operator
    SparseMatrix<std::complex<double>> vectorOp = MM + shortTime * LL;
    //SparseMatrix<std::complex<double>> vectorOp = id + shortTime * InvM*LL;

    // Check the Delaunay condition. If the mesh is Delaunay, then vectorOp is SPD, and we can use a
    // PositiveDefiniteSolver. Otherwise, we must use a SquareSolver
    geometry_bez_surf->requireEdgeCotanWeights();
    bool isDelaunay = false;// true;
    for (Edge e : mesh_bez_surf->edges()) {
        if (geometry_bez_surf->edgeCotanWeights[e] < -1e-6) {
            isDelaunay = false;
            break;
        }
    }
    geometry_bez_surf->unrequireEdgeCotanWeights();

    std::unique_ptr<LinearSolver<std::complex<double>>> vectorHeatSolver;
    if (isDelaunay) {
        vectorHeatSolver.reset(new PositiveDefiniteSolver<std::complex<double>>(vectorOp));
    }
    else {
        vectorHeatSolver.reset(new SquareSolver<std::complex<double>>(vectorOp)); // not necessarily SPD without Delaunay
    }

    //VectorHeatMethodSolver vhs(*geometry_bez_surf);
    geometry_bez_surf->requireVertexTangentBasis();
    const auto &vb = geometry_bez_surf->vertexTangentBasis;
    VertexData<Vector3> vbx(*mesh_bez_surf);
    VertexData<Vector3> vby(*mesh_bez_surf);
    VertexData<Eigen::dcomplex> q_vtxb(*mesh_bez_surf);
    VertexData<Eigen::dcomplex> q_cpx(*mesh_bez_surf);
    VertexData<double> q_mag(*mesh_bez_surf);
    VertexData<double> M(*mesh_bez_surf);
    for (auto v : mesh_bez_surf->vertices()) {
        vbx[v] = vb[v][0];
        vby[v] = vb[v][1];
        auto bs = Ns[v].buildTangentBasis();
        const auto& b1 = bs[0];
        const auto& b2 = bs[1];
        auto qv = Q[v][0]*bs[0] + Q[v][1]*bs[1];
        Vector2 q_cpx_v2 = { dot(qv, vb[v][0]), dot(qv, vb[v][1]) };
        q_vtxb[v] = q_cpx_v2;
        q_cpx_v2 = q_cpx_v2.pow(2);
        q_cpx[v] = q_cpx_v2;
        M[v] = MM.coeff(v.getIndex(), v.getIndex()).real();
        //q_cpx[v] *= MM.coeffRef(v.getIndex(), v.getIndex());
    }
    psMesh_bez_surf->addVertexVectorQuantity("vtx basis (x)", vbx);
    psMesh_bez_surf->addVertexVectorQuantity("vtx basis (y)", vby);
    psMesh_bez_surf->addVertexTangentVectorQuantity("Q (vtx basis)", q_vtxb, vbx, vby, 2);
    //psMesh_bez_surf->addVertexTangentVectorQuantity("LQ", InvM*(LL*q_cpx.toVector()), vbx, vby, 2);
    psMesh_bez_surf->addVertexScalarQuantity("Mass", M, polyscope::DataType::MAGNITUDE);

    int num_iter = num_smoothing_iters;
    while (num_iter-- > 0) {
        q_cpx = VertexData<std::complex<double>>(*mesh_bez_surf, vectorHeatSolver->solve(MM*q_cpx.toVector()));
    }

    //// Print
    //for (auto v : mesh_bez_surf->vertices()) {
    //    if (us[v] < 1E-6 && vs[v] < 1E-6) {
    //        Vector<Eigen::dcomplex> probe(mesh_bez_surf->nVertices());
    //        probe.setZero();
    //        probe(v.getIndex()) = 1.0;
    //        Vector<Eigen::dcomplex> res = InvM * (LL * probe);
    //        std::cerr << v.getIndex() << res(v.getIndex()) << std::endl;
    //        for (auto vvi : v.adjacentVertices()) {
    //            std::cerr << vvi.getIndex() << res(vvi.getIndex()) << std::endl;
    //        }

    //        //std::cerr << "Scaled angles: ";
    //        //for (auto sa : geometry_bez_surf->cornerScaledAngles[v.getIndex()]) {
    //        //    std::cerr << sa << ", ";
    //        //}
    //        //std::cerr << std::endl;
    //    }
    //}
    //std::cerr << "usesImplicitTwin:" << mesh_bez_surf->usesImplicitTwin() << std::endl;
    

    //VertexData<Vector2> q_cpx(*mesh_bez_surf);
    for (auto v : mesh_bez_surf->vertices()) { // Computing Q, S and II
        q_cpx[v] = Vector2::fromComplex(q_cpx[v]).pow(1.0/2.0);

        auto q_vec_1_3D = q_cpx[v].real() * vbx[v] + q_cpx[v].imag() * vby[v];
        auto q_vec_2_3D = -q_cpx[v].imag() * vbx[v] + q_cpx[v].real() * vby[v];
        q_mag[v] = q_vec_1_3D.norm();
        q_vec_1_3D = q_vec_1_3D.normalize();
        q_vec_2_3D = q_vec_2_3D.normalize();
        double H = mean[v];
        auto S_eval1 = H + q_mag[v];
        auto S_eval2 = H - q_mag[v];
        const auto& du = dus[v];
        const auto& dv = dvs[v];
        auto q_vec_1_uv = inSystem(q_vec_1_3D, du, dv);
        auto q_vec_2_uv = inSystem(q_vec_2_3D, du, dv);
        Eigen::Matrix2d EV, VV;
        EV(0, 0) = q_vec_1_uv[0]; EV(0, 1) = q_vec_2_uv[0];
        EV(1, 0) = q_vec_1_uv[1]; EV(1, 1) = q_vec_2_uv[1];
        VV = Eigen::Vector2d(S_eval1, S_eval2).asDiagonal();
        Eigen::Matrix2d SS = EV * (VV * EV.inverse());
        Eigen::Matrix2d I, II;
        double E = dot(du, du);
        double F = dot(du, dv);
        double G = dot(dv, dv);
        I(0, 0) = E; I(0, 1) = F;
        I(1, 0) = F; I(1, 1) = G;
        II = I * SS;
        double L = II(0, 0);
        double M = II(0, 1);
        double N = II(1, 1);

        auto duu_uv = duus[v].removeComponent(Ns[v]);
        duu_uv += L * Ns[v];
        auto duv_uv = duvs[v].removeComponent(Ns[v]);
        duv_uv += M * Ns[v];
        auto dvv_uv = dvvs[v].removeComponent(Ns[v]);
        dvv_uv += N * Ns[v];

        duus_smoothed[v] = duu_uv;
        duvs_smoothed[v] = duv_uv;
        dvvs_smoothed[v] = dvv_uv;


        //q_cpx[v] /= MM.coeffRef(v.getIndex(), v.getIndex());
    }
    psMesh_bez_surf->addVertexTangentVectorQuantity("Q (smoothed)", q_cpx, vbx, vby, 2);
    psMesh_bez_surf->addVertexScalarQuantity("Q (smoothed) magnitude", q_mag, polyscope::DataType::MAGNITUDE);
    

}

void smoothMeanCurvature()
{
    geometry_bez_surf->requireEdgeLengths();
    geometry_bez_surf->requireVertexGalerkinMassMatrix();
    geometry_bez_surf->requireCotanLaplacian();
    const auto &MM = geometry_bez_surf->vertexGalerkinMassMatrix;
    const auto &LL = geometry_bez_surf->cotanLaplacian;

    double tCoef = 1.0;
    double meanEdgeLength = 0.;
    for (Edge e : mesh_bez_surf->edges()) {
        meanEdgeLength += geometry_bez_surf->edgeLengths[e];
    }
    meanEdgeLength /= mesh_bez_surf->nEdges();
    double shortTime = tCoef* meanEdgeLength* meanEdgeLength;

    SparseMatrix<double> heatOp = MM + shortTime * LL;
    std::unique_ptr<PositiveDefiniteSolver<double>> scalarHeatSolver;
    scalarHeatSolver.reset(new PositiveDefiniteSolver<double>(heatOp));

    mean_smoothed = VertexData<double>(*mesh_bez_surf);
    auto mean_smoothed_eig = scalarHeatSolver->solve(MM*mean.toVector());
    for (auto v : mesh_bez_surf->vertices()) {
        mean_smoothed[v] = mean_smoothed_eig(v.getIndex());
    }
    psMesh_bez_surf->addVertexScalarQuantity("Mean curvature (smoothed)", mean_smoothed, polyscope::DataType::SYMMETRIC);
}

void fitSmoothedDerivatives()
{
    std::vector<Eigen::Vector3d> dse(duus_smoothed.size() * 3);
    std::vector<Eigen::Vector2d> uvs(duus_smoothed.size() * 3);
    for (size_t i = 0; i < duus_smoothed.size(); ++i) {
        dse[3 * i]     = { duus_smoothed[i][0], duus_smoothed[i][1], duus_smoothed[i][2] };
        dse[3 * i + 1] = { duvs_smoothed[i][0], duvs_smoothed[i][1], duvs_smoothed[i][2] };
        dse[3 * i + 2] = { dvvs_smoothed[i][0], dvvs_smoothed[i][1], dvvs_smoothed[i][2] };

        uvs[3 * i]     = { us[i], vs[i] };
        uvs[3 * i + 1] = { us[i], vs[i] };
        uvs[3 * i + 2] = { us[i], vs[i] };
    }
    std::vector<Eigen::Vector3d> new_cps;
    fitBeziertoDerivatives(dse, uvs, surf.basisU().degree(), surf.basisV().degree(), BSpline::toEigen(surf.controlPoints()), new_cps, 2, w_ders, w_oldcps);

    //surf.controlPoints() = BSpline::fromEigen(new_cps);
    //size_t resolution = 100;
    //BSpline::writeToMesh(surf, resolution);

    //std::tie(mesh_bez_surf, geometry_bez_surf) = readManifoldSurfaceMesh("bezier.obj");
    //psMesh_bez_cnet->updateVertexPositions(new_cps);

    psCurvNet_bez_smoothed_cnet->updateNodePositions(new_cps);
    polyscope::refresh();
}



// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {

  ImGui::SliderFloat("w_ders", &w_ders, 0.00000, 1.0, "%.5f");
  ImGui::SliderFloat("w_oldcps", &w_oldcps, 0.00000, 1.0, "%.5f");
  ImGui::SliderInt("# smoothing iters.", &num_smoothing_iters, 1, 100);

  if (ImGui::Button("Load Bezier")) {
      doWork();
  }

  if (ImGui::Button("Smooth derivatives")) {
      smoothDerivatives();
  }

  if (ImGui::Button("Fit smoothed derivatives")) {
      fitSmoothedDerivatives();
  }

  if (ImGui::Button("Smooth Weingartens")) {
      smoothWeingartens();
  }

  if (ImGui::Button("Smooth Qs")) {
      smoothQs();
  }

  if (ImGui::Button("Smooth mean curvature")) {
      smoothMeanCurvature();
  }

  if (ImGui::Button("Refresh normals")) {
      refreshNormals();
  }


}

int main(int argc, char **argv) {

  // Configure the argument parser
  args::ArgumentParser parser("geometry-central & Polyscope example project");
  args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help &h) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  //// Make sure a mesh name was given
  //if (!inputFilename) {
  //  std::cerr << "Please specify a mesh file as argument" << std::endl;
  //  return EXIT_FAILURE;
  //}

  // Initialize polyscope
  polyscope::init();
  polyscope::view::setUpDir(polyscope::UpDir::ZUp);
  polyscope::view::setFrontDir(polyscope::FrontDir::YFront);
  polyscope::options::ssaaFactor = 3;

  polyscope::loadStaticMaterial("isophotes", "isophotes.png");
  


  // Set the callback function
  polyscope::state::userCallback = myCallback;

  //// Load mesh
  //std::tie(mesh, geometry) = readManifoldSurfaceMesh(args::get(inputFilename));

  //// Register the mesh with polyscope
  //psMesh = polyscope::registerSurfaceMesh(
  //    polyscope::guessNiceNameFromPath(args::get(inputFilename)),
  //    geometry->inputVertexPositions, mesh->getFaceVertexList());

  //psMesh->setAllPermutations(polyscopePermutations(*mesh));

  // Set vertex tangent spaces
  //geometry->requireVertexTangentBasis();
  //VertexData<Vector3> vBasisX(*mesh), vBasisY(*mesh);
  //for (Vertex v : mesh->vertices()) {
  //  vBasisX[v] = geometry->vertexTangentBasis[v][0];
  //  vBasisY[v] = geometry->vertexTangentBasis[v][1];
  //}
  //psMesh->setVertexTangentBasisX(vBasisX);

  // auto vField =
  //     geometrycentral::surface::computeSmoothestVertexDirectionField(*geometry);
  // psMesh->addVertexTangentVectorQuantity("VF", vField, vBasisX, vBasisY);

//// Generate a guiding field
//  VertexData<Vector2> guideField =
//    geometrycentral::surface::computeSmoothestVertexDirectionField(*geometry, 2);
//
//  // Compute the stripe pattern
//  double constantFreq = 80.;
//  VertexData<double> frequencies(*mesh, constantFreq);
//  CornerData<double> periodicFunc;
//  FaceData<int> zeroIndices;
//  FaceData<int> branchIndices;
//  std::tie(periodicFunc, zeroIndices, branchIndices) =
//    computeStripePattern(*geometry, frequencies, guideField);
//
//  // Extract isolines
//  std::vector<Vector3> isolineVerts;
//  std::vector<std::array<size_t, 2>> isolineEdges;
//  std::tie(isolineVerts, isolineEdges) = extractPolylinesFromStripePattern(
//    *geometry, periodicFunc, zeroIndices, branchIndices, guideField, false);
//
//  geometry->requireVertexNormals();
//  psMesh->addVertexVectorQuantity("Vertex Normals", geometry->vertexNormals);
//
//  //for (auto v : mesh->vertices()) {
//  //    isolineVerts.push_back(geometry->vertexPositions[v]);
//  //}
//  //for (auto e : mesh->edges()) {
//  //    isolineEdges.push_back({e.firstVertex().getIndex(), e.secondVertex().getIndex() });
//  //}
//  polyscope::registerCurveNetwork("Stripe pattern", isolineVerts, isolineEdges);

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
