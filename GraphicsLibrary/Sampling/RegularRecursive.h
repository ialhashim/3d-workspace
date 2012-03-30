#include "GraphicsLibrary/Mesh/QSegMesh.h"
#include "MathLibrary/Bounding/BoundingBox.h"
#include "GraphicsLibrary/SpacePartition/Octree.h"
#include "GraphicsLibrary/SpacePartition/kdtree.h"

class RegularRecursiveSampling{

private:
	QSurfaceMesh * mesh;
	Octree * octree;
	KDTree * kd_tree;

	struct RRParam{
		double offset;
		double minDiag;
	};

	void SubdivideAndSample(std::vector<Vec3d> & pnts, BoundingBox bb, RRParam & rrp, float curDiag)
	{	
		Vec3d startPt = bb.Center();

		float dist_upper_bound = curDiag + rrp.offset;
		double dist = dist_upper_bound;

		// Compute mesh point nearest to bb center	
		IndexSet closestTris = octree->intersectPoint(startPt);
		
		Vec3d closestPt(0,0,0);
		double minDist = DBL_MAX;

		foreach(int idx, closestTris)
		{
			Vec3d q (mesh->closestPointFace(Surface_mesh::Face(idx), startPt));
			double curDist = (startPt - q).norm();

			if(curDist < minDist){
				closestPt = q;
				minDist = curDist;
			}
		}

		// = rrp.mesh.GetClosestPointDistance (startPt, dist_upper_bound, dist);

		curDiag /= 2;

		if(dist <= dist_upper_bound) 
		{
			// store points only for the last level of recursion (?)
			// Base case
			if(curDiag / 3 < rrp.minDiag) 
			{
				if(rrp.offset == 0)
				{
					if(!kd_tree->has(closestPt[0], closestPt[1], closestPt[2]) && closestPt.sqrnorm())
					{
						kd_tree->insert(&closestPt[0], 0);
						pnts.push_back(closestPt);
					}
				}
				else 
				{
					// points below the offset threshold cannot be displaced 
					// at the right offset distance, we can only make points nearer.
					if(dist > rrp.offset) 
					{
						Vec3d delta = startPt - closestPt;

						Vec3d q(closestPt + delta * (rrp.offset / dist));

						if(!kd_tree->has(q[0], q[1], q[2]) &&q.sqrnorm())
						{
							pnts.push_back( q );
							kd_tree->insert(&q[0], 0);
						}
					}
				}
			}

			if(curDiag < rrp.minDiag) return;

			Vec3d hs = (bb.vmax - bb.vmin) / 2;

			// Visit children
			for(int i = 0; i < 2; i++){
				for(int j = 0; j < 2; j++){
					for(int k = 0; k < 2; k++)
					{
						BoundingBox newBB( 
							Vec3d( bb.vmin[0]+i*hs[0], bb.vmin[1]+j*hs[1], bb.vmin[2]+k*hs[2]), 
							Vec3d(startPt[0]+i*hs[0],startPt[1]+j*hs[1],startPt[2]+k*hs[2]) );

						SubdivideAndSample(pnts, newBB, rrp, curDiag);
					}
				}
			}
		}
	} 

public:

	RegularRecursiveSampling(){}

	std::vector<Vec3d> sample(QSurfaceMesh * m, float minDiag, double offset = 0.0)
	{
		std::vector<Vec3d> pnts;

		this->mesh = m;
	
		// Collect mesh triangles, build octree
		std::list<BaseTriangle*> tris;
		for(int fi = 0; fi < m->n_faces(); fi++){
			std::vector<Vec3d> v = m->facePoints(Surface_mesh::Face(fi));
			tris.push_back(new Triangle(v[0],v[1],v[2], fi));
		}
		octree = new Octree;
		octree->initBuild(tris, 30);

		kd_tree = new KDTree;

		BoundingBox bb (m->bbmin, m->bbmax);
		bb.Offset(offset * 2.0);

		RRParam rrp;
		rrp.offset = offset;
		rrp.minDiag = minDiag;

		SubdivideAndSample( pnts, bb, rrp, bb.Diag() );

		return pnts;
	}

};
