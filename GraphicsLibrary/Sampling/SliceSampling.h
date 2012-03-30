#pragma once

#include "GraphicsLibrary/Mesh/QSegMesh.h"
#include "GraphicsLibrary/Basic/Plane.h"

class SliceSampling{
public:
	static std::vector<Vec3d> sample(QSurfaceMesh * m, int segments, double dt)
	{
		std::vector<Vec3d> pnts;

		Vec3d ZVEC(0,0,1);
		double deltaZ = (m->bbmax.z() - m->bbmin.z()) / segments;
		Vec3d delta(0,0, deltaZ);
		Vec3d center(0,0, m->bbmin.z() - 2 * deltaZ);

		// Compute contours
		for(int i = 0; i <= segments + 2; i++)
		{
			Plane slicePlane(ZVEC, center);

			std::vector<Vec3d> sliceSamples;

			for(int f = 0; f < m->n_faces(); f++){
				std::vector<Point> fv = m->facePoints(Surface_mesh::Face(f));

				Vec3d p1(0,0,0), p2(0,0,0);

				int RESULT = slicePlane.ContourFacet(fv[0], fv[1], fv[2], p1, p2);

				if(RESULT == INTERSECT || RESULT == ENDPOINT_INTERSECT)
				{
					if(p1.sqrnorm()) sliceSamples.push_back(p1);
					if(p2.sqrnorm()) sliceSamples.push_back(p2);
				}
			}

			// Add to result set
			foreach(Point p, sliceSamples)
				pnts.push_back(p);

			center += delta;
		}

		return pnts;
	}

	static std::vector<Vec3d> sample2(QSurfaceMesh * m, int segments, double dt)
	{
		std::vector<Vec3d> pnts;

		double width = m->bbmax.x() - m->bbmin.x();
		double length = m->bbmax.y() - m->bbmin.y();
		double height = m->bbmax.z() - m->bbmin.z();

		Vec3d ZVEC(0,0,1);
		double deltaZ = height / segments;

		Vec3d corner(m->bbmin.x(), m->bbmin.y(), 0);
		Vec3d bottom(0, 0, m->bbmin.z());

		int dx = width / dt;
		int dy = length / dt;

		double r2 = dt * 0.1;

		std::vector< std::vector<Vec3d> > grid;

		for(int y = 0; y <= dy; y++){
			grid.push_back(std::vector<Vec3d>());
			for(int x = 0; x <= dx; x++)
				grid[y].push_back(corner + Vec3d((double(x) / dx) * width, (double(y) / dy) * length, bottom.z()));
		}

		// Compute contours
		for(int i = 0; i <= segments; i++)
		{
			Plane slicePlane(ZVEC, bottom + Vec3d(0,0,(double(i)/segments) * height));

			std::vector<Vec3d> sliceSamples;

			// Collect triangles intersecting planes
			std::set<int> tris;
			for(int f = 0; f < m->n_faces(); f++){
				std::vector<Point> fv = m->facePoints(Surface_mesh::Face(f));
				
				int vtest0 = 0, vtest1 = 0, vtest2 = 0;

				if(slicePlane.IsFront(fv[0])) vtest0 = 1;
				if(slicePlane.IsFront(fv[1])) vtest1 = 1;
				if(slicePlane.IsFront(fv[2])) vtest2 = 1;

				int test = vtest0 + vtest1 + vtest2;

				if(test < 3)
					tris.insert(f);
			}

			// Check grid points intersecting triangles
			foreach(int f, tris)
			{
				std::vector<Point> fv = m->facePoints(Surface_mesh::Face(f));

				for(int y = 0; y < grid.size(); y++)
				{
					for(int x = 0; x < grid.front().size(); x++)
					{
						Vec3d center = grid[y][x] + Vec3d(0,0,(double(i)/segments) * height);

						Vec3d p = ClosestPtPointTriangle2(center,  fv[0],  fv[1],  fv[2]);

						double dist = (p - center).norm();

						if (dist <= r2)
							sliceSamples.push_back(p);
					}
				}
			}

			// Add to result set
			foreach(Point p, sliceSamples)
				pnts.push_back(p);
		}

		return pnts;
	}
};
