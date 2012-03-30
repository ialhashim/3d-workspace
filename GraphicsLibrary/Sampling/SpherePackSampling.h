#pragma once

#include "Sampler.h"
#include "GraphicsLibrary/SpacePartition/kdtree.h"

class SpherePackSampling{

public:
	static std::vector<Vec3d> sample(QSurfaceMesh * m, int randomSampleCount, double r, int density = 1)
	{
		std::vector<Vec3d> samples, centers, rndSamples;

		// Centers of packed spheres
		centers = spheres(r, m->bbmin, m->bbmax, density);

		// Get a lot of random samples
		foreach(SamplePoint sp, Sampler(m).getSamples(randomSampleCount)){
			rndSamples.push_back(sp.pos);
		}

		// Add into a KD-tree
		KDTree tree;
		foreach(Point p, rndSamples) tree.insert(&p[0], 0);

		foreach(Point center, centers)
		{
			// Collect neighbors
			kdres * neighbors = tree.nearest_range(&center[0], r);
			res_node * node = neighbors->riter;

			if(node == NULL) continue;

			std::vector<Vec3d> group;

			while(node != NULL)	{
				double * pos = node->item->pos;
				group.push_back(Vec3d(pos[0], pos[1], pos[2]));
				node = node->next;
			}

			// Record center
			Vec3d centerGroup (0,0,0);
			foreach(Point p, group) centerGroup += p;
			centerGroup /= group.size();
			samples.push_back(centerGroup);
		}

		return samples;
	}

	/* Hexagonal close packing of spheres (HCP lattice) */
	static std::vector<Vec3d> spheres(double r, Vec3d bbmin, Vec3d bbmax, int density = 1)
	{
		std::vector<Vec3d> samples;

		std::vector< std::vector<Vec3d> > grid;

		double d = r * 2;

		double width = bbmax.x() - bbmin.x();
		double length = bbmax.y() - bbmin.y();
		double height = bbmax.z() - bbmin.z();

		Vec3d corner(bbmin.x()-r, bbmin.y()-r, bbmin.z());

		int dx = (width / r) + 1;
		int dy = (length / r) + 1;
		int dz = (height / r) + 1;

		Vec3d delta(r, sqrt(3.0) * r , 0);
		Vec3d center(0,0,0);
		int sign = 1;

		grid.push_back(std::vector<Vec3d>());

		for(int x = 0; x < dx; x++)
			grid[0].push_back(corner + Vec3d(x * d, 0, 0));

		for(int y = 0; y < dy; y++){
			std::vector<Vec3d> row = grid.back();
			
			for(int x = 0; x < row.size(); x++)	row[x] += delta;

			delta.x() = -delta.x();

			grid.push_back(row);
		}

		Vec3d omega(r, r, sqrt(6.0) * (2.0/3.0) * r);

		for(int z = 0; z < dz; z += density){
			for(int y = 0; y < dy; y += density)
				for(int x = 0; x < dx; x += density)
					samples.push_back(grid[y][x] + center);
				
			omega.x() *= -1;
			omega.y() *= -1;

			center += omega;
		}

		return samples;
	}

};
