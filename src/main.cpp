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

// Some algorithm parameters
float param1 = 42.0;

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

    Geometry::BSSurface surf(deg_u, deg_v, cps);

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

    //psMesh_bez_cnet->setAllPermutations(polyscopePermutations(*mesh_bez_cnet));


    VertexData<double> us(*mesh_bez_surf);
    VertexData<double> vs(*mesh_bez_surf);

    VertexData<Vector3> dus(*mesh_bez_surf);
    VertexData<Vector3> dvs(*mesh_bez_surf);
    VertexData<Vector3> duus(*mesh_bez_surf);
    VertexData<Vector3> duvs(*mesh_bez_surf);
    VertexData<Vector3> dvvs(*mesh_bez_surf);

    VertexData<double> mean(*mesh_bez_surf);
    VertexData<double> gauss(*mesh_bez_surf);
    VertexData<double> DJ(*mesh_bez_surf);
    VertexData<Vector2> dmin(*mesh_bez_surf);
    VertexData<Vector2> dmax(*mesh_bez_surf);
    VertexData<Vector2> Q(*mesh_bez_surf);

    VertexData<Vector3> Ns(*mesh_bez_surf);

    geometry_bez_surf->requireVertexNormals();
    geometry_bez_surf->requireVertexDualAreas();

    psMesh_bez_surf->vertexNormals.ensureHostBufferAllocated();

    for (size_t i = 0; i < resolution; ++i) {
        double u = (double)i / (double)(resolution - 1);
        for (size_t j = 0; j < resolution; ++j) {
            double v = (double)j / (double)(resolution - 1);
            
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

            Eigen::Matrix2d I;
            I(0, 0) = E; I(0, 1) = F;
            I(1, 0) = F; I(1, 1) = G;

            Eigen::Matrix2d II;
            II(0, 0) = L; II(0, 1) = M;
            II(1, 0) = M; II(1, 1) = N;

            Eigen::Matrix2d S = II * I.inverse();

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

            Ns[v_idx] = { -n_[0], -n_[1], -n_[2] };
            DJ[v_idx] = n.norm();

            Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es_S(S);
            Eigen::Vector2d e1_uv = es_S.eigenvectors().col(0);
            Eigen::Vector2d e2_uv = es_S.eigenvectors().col(1);
            //auto e1_xyz = du * e1_uv[0] + dv * e1_uv[1];
            //auto e2_xyz = du * e1_uv[0] + dv * e1_uv[1];
            //dmin[v_idx] = { e1_xyz[0], e1_xyz[1], e1_xyz[2] };
            //dmax[v_idx] = { e2_xyz[0], e2_xyz[1], e2_xyz[2] };
            dmin[v_idx] = { e1_uv[0], e1_uv[1] };
            dmax[v_idx] = { e2_uv[0], e2_uv[1] };
            Eigen::Matrix2d Hopf = S - Eigen::Matrix2d::Identity() * H;
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es_Q(Hopf);
            e1_uv = std::abs(es_Q.eigenvalues()[0])*es_Q.eigenvectors().col(0);
            e2_uv = std::abs(es_Q.eigenvalues()[1])*es_Q.eigenvectors().col(1);
            Q[v_idx] = { e1_uv[0], e1_uv[1] };

            geometry_bez_surf->vertexNormals[v_idx] = Ns[v_idx];

        }
    }

    VertexData<Vector3> duus_smoothed(*mesh_bez_surf);
    //VectorHeatMethodSolver vhs(*geometry_bez_surf);
    //geometry_bez_surf->requireCotanLaplacian();
    //geometry_bez_surf->requireVertexLumpedMassMatrix();
    //const auto& LL = geometry_bez_surf->cotanLaplacian;
    //const auto& MM = geometry_bez_surf->vertexLumpedMassMatrix;
    //double dt = 1.0;
    //SparseMatrix<double> heat = MM + dt*LL;
    //PositiveDefiniteSolver solv(heat);
    //for (size_t cc = 0; cc < 3; ++cc) {
    //    VertexData<double> rhs(*mesh_bez_surf);
    //    for (auto v : mesh_bez_surf->vertices()) {
    //        rhs[v] = duus[v][cc];
    //    }
    //    //auto MM = geometry_bez_surf->vertexLumpedMassMatrix;
    //    //
    //    //VertexData<double> x_pw(*mesh_bez_surf, solv.solve(x.toVector()));
    //    //VertexData<double> x(*mesh_bez_surf, solv.solve(MM*(rhs.toVector())));
    //    auto x = rhs;
    //    size_t num_iter = 10;
    //    while (num_iter-- > 0) {
    //        x = vhs.scalarDiffuse(VertexData<double>(*mesh_bez_surf,MM*x.toVector()));
    //    }

    //    //VertexData<double> x = vhs.scalarDiffuse());
    //    for (auto v : mesh_bez_surf->vertices()) {
    //        duus_smoothed[v][cc] = x[v];
    //    }
    //}


    psMesh_bez_surf->addVertexScalarQuantity("u", us, polyscope::DataType::MAGNITUDE);
    psMesh_bez_surf->addVertexScalarQuantity("v", vs, polyscope::DataType::MAGNITUDE);
    psMesh_bez_surf->addVertexVectorQuantity("Du", dus);
    psMesh_bez_surf->addVertexVectorQuantity("Dv", dvs);
    psMesh_bez_surf->addVertexVectorQuantity("Duu", duus);
    psMesh_bez_surf->addVertexVectorQuantity("Duu (smoothed)", duus_smoothed);
    psMesh_bez_surf->addVertexVectorQuantity("Duv", duvs);
    psMesh_bez_surf->addVertexVectorQuantity("Dvv", dvvs);
    psMesh_bez_surf->addVertexVectorQuantity("Vertex normals (mesh)", geometry_bez_surf->vertexNormals);
    psMesh_bez_surf->addVertexVectorQuantity("Vertex normals (exact)", Ns);
    psMesh_bez_surf->addVertexTangentVectorQuantity("dmin", dmin, dus, dvs);
    psMesh_bez_surf->addVertexTangentVectorQuantity("dmax", dmax, dus, dvs);
    psMesh_bez_surf->addVertexTangentVectorQuantity("Q", Q, dus, dvs, 2);
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
    psMesh_bez_surf->addVertexScalarQuantity("DJ",
        DJ,
        polyscope::DataType::MAGNITUDE);
    //psMesh->addEdgeScalarQuantity("edgelength",
    //    geometry->edgeLengths,
    //    polyscope::DataType::MAGNITUDE);

    psMesh_bez_surf->setShadeStyle(polyscope::MeshShadeStyle::Smooth);
    psMesh_bez_surf->refresh();
    polyscope::view::resetCameraToHomeView();
}

void refreshNormals() {

    for (auto v : mesh_bez_surf->vertices()) {
        size_t v_idx = v.getIndex();
        auto Ns = static_cast<polyscope::SurfaceVertexVectorQuantity*>(psMesh_bez_surf->getQuantity("Vertex normals (exact)"));
        psMesh_bez_surf->vertexNormals.data[v_idx] = Ns->vectors.data[v_idx] ;

    }
    psMesh_bez_surf->vertexNormals.markHostBufferUpdated();
    //psMesh_bez_surf->refresh();
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {

  if (ImGui::Button("Load Bezier")) {
      doWork();
  }

  if (ImGui::Button("Refresh normals")) {
      refreshNormals();
  }

  //ImGui::SliderFloat("param", &param1, 0., 100.);
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
