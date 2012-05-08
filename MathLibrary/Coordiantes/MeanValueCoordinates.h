#pragma once

#include "GraphicsLibrary/Mesh/QSurfaceMesh.h"

class MeanValueCooridnates{

public:
	static std::vector<double> weights(const Point& x, QSurfaceMesh * mesh)
	{
		int npts = mesh->n_vertices();

		// Initializing weights
		std::vector<double> weights(npts);

		// arrays storing point-to-vertex vectors and distances
		std::vector<double> dist(npts);
		std::vector<Vec3d> uVec(npts);
		static const double eps = 0.00000001;

		Surface_mesh::Vertex_property<Point> mesh_points = mesh->vertex_property<Point>("v:point");
		Surface_mesh::Vertex_iterator vit, vend = mesh->vertices_end();

		for (vit = mesh->vertices_begin(); vit != vend; ++vit)
		{
			int pid = Surface_mesh::Vertex(vit).idx();

			// point-to-vertex vector
			uVec[pid] = mesh_points[vit] - x;

			// distance
			dist[pid] = uVec[pid].norm();

			// handle special case when the point is really close to a vertex
			if (dist[pid] < eps){
				weights[pid] = 1.0;
				return weights;
			}

			// project onto unit sphere (normalize)
			uVec[pid] /= dist[pid];
		}

		if(!mesh->triangles.size()) mesh->fillTrianglesList();

		// Now loop over all triangle to compute weights
		for(int fi = 0; fi < mesh->triangles.size() / 3; fi++)
		{
			// vertex id
			int pid0 = mesh->triangles[3*fi + 0];
			int pid1 = mesh->triangles[3*fi + 1];
			int pid2 = mesh->triangles[3*fi + 2];

			// unit vectors
			Vec3d u0 (uVec[pid0]);
			Vec3d u1 (uVec[pid1]);
			Vec3d u2 (uVec[pid2]);

			// edge lengths
			double l0 = (u1 - u2).norm();
			double l1 = (u2 - u0).norm();
			double l2 = (u0 - u1).norm();

			// angles
			double theta0 = 2.0*asin(l0/2.0);
			double theta1 = 2.0*asin(l1/2.0);
			double theta2 = 2.0*asin(l2/2.0);
			double halfSum = (theta0 + theta1 + theta2) / 2.0;

			// special case when the point lies on the triangle
			if (M_PI - halfSum < eps)
			{
				weights.clear();
				weights.resize(npts, 0.0); // clear all

				weights[pid0] = sin(theta0) * dist[pid1] * dist[pid2];
				weights[pid1] = sin(theta1) * dist[pid2] * dist[pid0];
				weights[pid2] = sin(theta2) * dist[pid0] * dist[pid1];

				double sumWeight = weights[pid0] + weights[pid1] + weights[pid2];

				weights[pid0] /= sumWeight;
				weights[pid1] /= sumWeight;
				weights[pid2] /= sumWeight;

				return weights;
			}

			// coefficient
			double sinHalfSum = sin(halfSum);
			double sinHalfSumSubTheta0 = sin(halfSum-theta0);
			double sinHalfSumSubTheta1 = sin(halfSum-theta1);
			double sinHalfSumSubTheta2 = sin(halfSum-theta2);
			double sinTheta0 = sin(theta0), sinTheta1 = sin(theta1), sinTheta2 = sin(theta2);

			double c0 = 2 * sinHalfSum * sinHalfSumSubTheta0 / sinTheta1 / sinTheta2 - 1;
			double c1 = 2 * sinHalfSum * sinHalfSumSubTheta1 / sinTheta2 / sinTheta0 - 1;
			double c2 = 2 * sinHalfSum * sinHalfSumSubTheta2 / sinTheta0 / sinTheta1 - 1;

			if (fabs(c0) > 1) c0 = c0 > 0 ? 1 : -1;
			if (fabs(c1) > 1) c1 = c1 > 0 ? 1 : -1;
			if (fabs(c2) > 1) c2 = c2 > 0 ? 1 : -1;

			// sign
			double det = Determinant3x3(u0, u1, u2);

			// skip when less than eps
			if (abs(det) < eps){
				fi++; continue;
			}

			double detSign = det > 0 ? 1 : -1;
			double sign0 = detSign * sqrt(1 - c0*c0);
			double sign1 = detSign * sqrt(1 - c1*c1);
			double sign2 = detSign * sqrt(1 - c2*c2);

			// if 'x' lies on the plane of current triangle but outside it, ignore the current triangle
			if (abs(sign0) < eps || abs(sign1) < eps || abs(sign2) < eps)
			{
				fi++; continue;
			}

			// weight 
			weights[pid0] += (theta0-c1*theta2-c2*theta1) / (dist[pid0]*sinTheta1*sign2);
			weights[pid1] += (theta1-c2*theta0-c0*theta2) / (dist[pid1]*sinTheta2*sign0);
			weights[pid2] += (theta2-c0*theta1-c1*theta0) / (dist[pid2]*sinTheta0*sign1);
		}

		// normalize weight
		double sumWeight = 0.0;
		for (int pid=0; pid < npts; ++pid)	sumWeight += weights[pid];
		if(!sumWeight) printf("WARNING: zero weights.\n");
		for (int pid=0; pid < npts; ++pid)	weights[pid] /= sumWeight;

		return weights;
	}

	static Vec3d point( std::vector<double> weight, QSurfaceMesh * mesh )
	{
		Surface_mesh::Vertex_property<Point> mesh_points = mesh->vertex_property<Point>("v:point");
		Surface_mesh::Vertex_iterator vit, vend = mesh->vertices_end();

		Vec3d p(0,0,0);

		for (vit = mesh->vertices_begin(); vit != vend; ++vit){
			p += mesh_points[vit] * weight[Surface_mesh::Vertex(vit).idx()];
		}

		return p;
	}

	static inline double Determinant3x3(const Vec3d& c1, const Vec3d& c2, const Vec3d& c3)
	{
		return c1[0] * c2[1] * c3[2] + c2[0] * c3[1] * c1[2] + c3[0] * c1[1] * c2[2] -
			c1[0] * c3[1] * c2[2] - c2[0] * c1[1] * c3[2] - c3[0] * c2[1] * c1[2];
	}
};
