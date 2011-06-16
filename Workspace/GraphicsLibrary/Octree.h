// Based on Radakan engine (ported from Java) (GNU General Public License v3)
// With major modifications by Ibraheem
// http://radakan.googlecode.com/svn/trunk/simplephysics/src/com/jme/bounding/Octree.java

#pragma once

#include <cmath>

#include "BoundingBox.h"
#include "Triangle.h"

typedef std::set<int> IndexSet;
typedef IndexSet::iterator IndexSetIter;

#include <stack>
using namespace std;

class Octree
{
private:
	Vector<Octree> children;
	Vector<BaseTriangle*> triangleData;

public:
	BoundingBox boundingBox;
	int trianglePerNode;

	Octree(){trianglePerNode = -1;}
//	Octree(const Vector<int> & trisIndex, Mesh * mesh, int triPerNode);
        Octree(const StdList<BaseTriangle*>& tris, int triPerNode);
        Octree( int triPerNode, const BoundingBox& bb, const Vector<BaseTriangle*>& tris );

	void init(int triPerNode);
	void initBuild(const StdList<BaseTriangle*>& tris, int triPerNode );

	void newNode( int depth, double x, double y, double z );
	void build(int depth = 0);

        vector<BaseTriangle*> getIntersectingTris(const Vec& v0, const Vec& v1, const Vec& v2, bool showIt=false);

	bool intersectHit(IndexSet& tris);

	IndexSet intersectPoint(const Vec& point);
	void intersectRecursivePoint(const Vec& point, IndexSet& tris);

	IndexSet intersectRay(const Ray& ray);
	void intersectRecursiveRay(const Ray& ray, IndexSet& tris);

	IndexSet intersectSphere(const Vec& sphere_center, double radius);
	void intersectRecursiveSphere(const Vec& sphere_center, double radius, IndexSet& tris);

	IndexSet intersectRaySphere(const Ray& ray, const Vec& sphere_center, double radius);
	void intersectRayBoth(const Ray& ray, IndexSet & tris);

//	BaseTriangle* findClosestTri(const Ray & ray, IndexSet & tris, Mesh * mesh, HitResult & hitRes);

        /* Perform intersection tests  */
	bool testIntersectHit(const Ray& ray, HitResult & hitRes);
	void testIntersectRayBoth(const Ray& ray, HitResult & hitRes);

	Vector<BaseTriangle*> getTriangleData();

	void draw(double r, double g, double b);
};

