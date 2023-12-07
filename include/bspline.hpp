#pragma once

#include "geometry.hh"
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace BSpline{
	using namespace Geometry;
    
    size_t toIdx(size_t row, size_t col, size_t m, size_t n) 
    {
        return m*row + col;
    }

    std::tuple<size_t, size_t> fromIdx(size_t idx, size_t m, size_t n) 
    {
        return{ idx / m, idx % m };
    }

    Eigen::Vector3d toEigen(const Vector3D &gv) 
    {
        return { gv[0], gv[1], gv[2] };
    }

    Vector3D fromEigen(const Eigen::Vector3d &ev) 
    {
        return { ev[0], ev[1], ev[2] };
    }

    std::vector<Eigen::Vector3d> toEigen(const PointVector &gvs) 
    {
        std::vector<Eigen::Vector3d> evs(gvs.size());
        for (size_t i = 0; i < evs.size(); ++i) {
            evs[i] = { gvs[i][0], gvs[i][1], gvs[i][2] };
        }
        return evs;
    }

    PointVector fromEigen(const std::vector<Eigen::Vector3d> &evs) 
    {
        PointVector gvs(evs.size());
        for (size_t i = 0; i < gvs.size(); ++i) {
            gvs[i] = { evs[i][0], evs[i][1], evs[i][2] };
        }
        return gvs;
    }

    void toEigen(const Vector3D &gv, Eigen::Vector3d &ev) 
    {
        ev = toEigen(gv);
    }

    void fromEigen(const Eigen::Vector3d &ev, Vector3D &gv) 
    {
        gv = fromEigen(ev);
    }

    PointVector zeroPoints(const PointVector &same_as) 
    {
        return PointVector(same_as.size(), Vector3D(0,0,0));
    }


	bool writeToMesh(const BSSurface &surf, size_t resolution) 
    {
        TriMesh trimesh;

        PointVector pv, tri;
        size_t n = surf.basisU().degree(), m = surf.basisV().degree();

        std::vector<double> coeff_u, coeff_v;
        for (size_t i = 0; i < resolution; ++i) {
            double u = (double)i / (double)(resolution - 1);
            for (size_t j = 0; j < resolution; ++j) {
                double v = (double)j / (double)(resolution - 1);
                pv.push_back(surf.eval(u, v));
            }
        }
        trimesh.setPoints(pv);


        for (size_t i = 0; i < resolution - 1; ++i) {
            for (size_t j = 0; j < resolution - 1; ++j) {
                trimesh.addTriangle(
                    i * resolution + j,
                    i * resolution + j + 1,
                    (i + 1) * resolution + j
                );
                trimesh.addTriangle(
                    (i + 1) * resolution + j,
                    i * resolution + j + 1,
                    (i + 1) * resolution + j + 1
                );
            }
        }

        trimesh.writeOBJ("bezier.obj");
        return true;
	}
}