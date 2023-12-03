#pragma once

#include "geometry.hh"

namespace BSpline{
	using namespace Geometry;
	bool writeToMesh(const BSSurface& surf, size_t resolution) {
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