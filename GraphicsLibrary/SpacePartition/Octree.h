// Based on Radakan engine (ported from Java) (GNU General Public License v3)
// With major modifications by Ibraheem
// http://radakan.googlecode.com/svn/trunk/simplephysics/src/com/jme/bounding/Octree.java

#pragma once

#include <cmath>

#include "MathLibrary/Bounding/BoundingBox.h"
#include "GraphicsLibrary/Basic/Triangle.h"

typedef std::set<int> IndexSet;
typedef IndexSet::iterator IndexSetIter;

#include <stack>
using namespace std;

class Octree
{
private:
	StdVector<Octree> children;
	StdVector<BaseTriangle*> triangleData;

public:
	BoundingBox boundingBox;
	int trianglePerNode;

	Octree(){trianglePerNode = -1; parent = NULL;}

    Octree(StdList<BaseTriangle*>& tris, int triPerNode);
    Octree( int triPerNode, const BoundingBox& bb, const StdVector<BaseTriangle*>& tris );

	void init(int triPerNode);
	void initBuild(StdList<BaseTriangle*>& tris, int triPerNode );

	void newNode( int depth, double x, double y, double z );
	void build(int depth = 0);

    StdVector<BaseTriangle*> getIntersectingTris(const Vec3d& v0, const Vec3d& v1, const Vec3d& v2, bool showIt=false);

	bool intersectHit(IndexSet& tris);

	IndexSet intersectPoint(const Vec3d& point);
	void intersectRecursivePoint(const Vec3d& point, IndexSet& tris);

	IndexSet intersectRay(const Ray& ray);
	void intersectRecursiveRay(const Ray& ray, IndexSet& tris);

	IndexSet intersectSphere(const Vec3d& sphere_center, double radius);
	void intersectRecursiveSphere(const Vec3d& sphere_center, double radius, IndexSet& tris);

	IndexSet intersectRaySphere(const Ray& ray, const Vec3d& sphere_center, double radius);
	void intersectRayBoth(const Ray& ray, IndexSet & tris);

//	BaseTriangle* findClosestTri(const Ray & ray, IndexSet & tris, Mesh * mesh, HitResult & hitRes);

        /* Perform intersection tests  */
	bool testIntersectHit(const Ray& ray, HitResult & hitRes);
	void testIntersectRayBoth(const Ray& ray, HitResult & hitRes);

	StdVector<BaseTriangle*> getTriangleData();

	Octree * parent;
	Octree * root();

	void draw(double r, double g, double b, double lineWidth = 1.0);

	// Debug
	std::vector<Octree *> selectedChildren;
};

