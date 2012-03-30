#include "Octree.h"
#include "Utility/SimpleDraw.h"

//Octree::Octree( const StdVector<int> & trisIndex, Mesh * mesh, int triPerNode )
//{
//	for(StdVector<int>::const_iterator index = trisIndex.begin(); index != trisIndex.end(); index++)
//	{
//		this->triangleData.push_back(mesh->f(*index));
//	}

//	init(triPerNode);
//}

Octree::Octree( StdList<BaseTriangle*>& tris, int triPerNode )
{
	for(StdList<BaseTriangle*>::const_iterator f = tris.begin(); f != tris.end(); f++)
	{
		this->triangleData.push_back(*f);
	}

	init(triPerNode);
}

Octree::Octree( int triPerNode, const BoundingBox& bb, const StdVector<BaseTriangle*>& tris )
{
	this->boundingBox = bb;
	this->trianglePerNode = triPerNode;
	this->triangleData = tris;
}

void Octree::init( int triPerNode )
{
	this->trianglePerNode = triPerNode;

	BoundingBox bb;
	bb.computeFromTris(triangleData);

	//Transform and scale to node's coordinates
	double largeSize = Max(bb.xExtent, Max(bb.yExtent, bb.zExtent));

	this->boundingBox = BoundingBox(bb.center, largeSize, largeSize, largeSize);

	parent = NULL;
}

void Octree::initBuild( StdList<BaseTriangle*>& tris, int triPerNode )
{
	// add triangles involved to "triangleData"
	this->triangleData = LIST_TO_VECTOR(tris);
	this->trianglePerNode = triPerNode;

	// Create a big box
	BoundingBox bb;
	bb.computeFromTris(triangleData);

	// Transform and scale to node's coordinates
	double largeSize = Max(bb.xExtent, Max(bb.yExtent, bb.zExtent));

	// Define our bounding box
	this->boundingBox = BoundingBox(bb.center, largeSize, largeSize, largeSize);

	// Build the tree
	this->build();

	// Connect children with parent
	std::stack<Octree*> childStack;
	childStack.push(this);
	while(!childStack.empty())
	{
		Octree * curr = childStack.top(); childStack.pop();

		for(int i = 0; i < (int) curr->children.size(); i++)
		{
			curr->children[i].parent = curr;

			childStack.push(&curr->children[i]);
		}
	}
}

void Octree::newNode( int depth, double x, double y, double z )
{
	double extent = boundingBox.xExtent / 2.0;

	Vec3d center;

	center.x() = boundingBox.center.x() + (extent * x);
	center.y() = boundingBox.center.y() + (extent * y);
	center.z() = boundingBox.center.z() + (extent * z);

	BoundingBox bb(center, extent, extent, extent);

	// Add child
	children.push_back(Octree());
	Octree * child = &children.back();

	child->boundingBox = bb;
	child->trianglePerNode = this->trianglePerNode;

	// Collect triangles inside child's bounding box
	for(StdVector<BaseTriangle*>::iterator it = this->triangleData.begin(); it != this->triangleData.end(); it++)
	{
		BaseTriangle* face = *it;

		if( bb.containsTriangle(face->vec(0), face->vec(1), face->vec(2)) )
		{
			child->triangleData.push_back(face);
		}
	}

	child->build(depth + 1); // build it
}

StdVector<BaseTriangle*> Octree::getIntersectingTris(const Vec3d& v0, const Vec3d& v1, const Vec3d& v2, bool showIt)
{
	if(this->triangleData.size() == 0 || this->children.size() == 0)
		return this->getTriangleData();

	StdVector<BaseTriangle*> res;

	for(int i = 0; i < (int) this->children.size(); i++)
	{
		if(children[i].boundingBox.containsTriangle(v0, v1, v2))
		{
			const vector <BaseTriangle*> &tris = children[i].getIntersectingTris(v0, v1, v2, showIt);

			for(int j = 0; j < (int) tris.size(); j++)
				res.push_back(tris[j]);
		}
	}

	return res;
}

void Octree::build( int depth /*= 0*/ )
{
	if ((int)triangleData.size() > this->trianglePerNode)
	{
		if(depth < 10)
		{
			// Subdivide to 8 nodes
			newNode(depth, -1, -1, -1);
			newNode(depth, 1, -1, -1);
			newNode(depth, -1, 1, -1);
			newNode(depth, 1, 1, -1);
			newNode(depth, -1, -1, 1);
			newNode(depth, 1, -1, 1);
			newNode(depth, -1, 1, 1);
			newNode(depth, 1, 1, 1);
		}
	}
}

void Octree::draw( double r, double g, double b, double lineWidth )
{
	SimpleDraw::DrawBox(boundingBox.center, boundingBox.xExtent, boundingBox.yExtent, boundingBox.zExtent,r,g,b, lineWidth);

	for (StdVector<Octree>::iterator child = children.begin();  child != children.end(); child++)
		child->draw(r,g,b, lineWidth);
}

IndexSet Octree::intersectPoint( const Vec3d& point )
{
	IndexSet tris;

	if (boundingBox.contains(point)) 
		intersectRecursivePoint(point, tris);

	return tris;
}

void Octree::intersectRecursivePoint( const Vec3d& point, IndexSet& tris )
{
	if (intersectHit(tris))
		return;

	for (StdVector<Octree>::iterator child = children.begin();  child != children.end(); child++)
	{
		if (child->boundingBox.contains(point))
			child->intersectRecursivePoint(point, tris);
	}
}

bool Octree::intersectHit( IndexSet& tris )
{
	if(this->children.size() > 0)
		return false;

	for(StdVector<BaseTriangle*>::iterator it = triangleData.begin(); it != triangleData.end(); it++)
	{
		BaseTriangle * face = *it;

		tris.insert( face->index );
	}

	root()->selectedChildren.push_back(this);

	return true;
}

IndexSet Octree::intersectRay( const Ray& ray )
{
	IndexSet tris;

	if (boundingBox.intersects(ray)) 
		intersectRecursiveRay(ray, tris);

	return tris;
}

void Octree::intersectRecursiveRay( const Ray& ray, IndexSet& tris )
{
	if (intersectHit(tris))
		return;

	for (StdVector<Octree>::iterator child = children.begin();  child != children.end(); child++)
	{
		if (child->boundingBox.intersects(ray))
			child->intersectRecursiveRay(ray, tris);
	}
}

IndexSet Octree::intersectSphere( const Vec3d& sphere_center, double radius )
{
	IndexSet tris;

	if (boundingBox.intersectsSphere(sphere_center, radius)) 
		intersectRecursiveSphere(sphere_center, radius, tris);

	return tris;
}

void Octree::intersectRecursiveSphere( const Vec3d& sphere_center, double radius, IndexSet& tris )
{
	if (intersectHit(tris))
		return;

	for (StdVector<Octree>::iterator child = children.begin();  child != children.end(); child++)
	{
		if (child->boundingBox.intersectsSphere(sphere_center, radius))
			child->intersectRecursiveSphere(sphere_center, radius, tris);
	}
}

IndexSet Octree::intersectRaySphere( const Ray& ray, const Vec3d& sphere_center, double radius )
{
	IndexSet tris;

	intersectRayBoth(ray, tris);

	// Grow sphere until we intersect ?
	//while(!boundingBox.intersectsSphere(sphere_center, radius))
	//	radius *= 1.25f;

	if(boundingBox.intersectsSphere(sphere_center, radius))
		intersectRecursiveSphere(sphere_center, radius, tris);

	return tris;
}

void Octree::intersectRayBoth( const Ray& ray, IndexSet & tris )
{
	stack<Octree*> trees;

	Ray inverseRay(ray.inverse());

	if (boundingBox.intersects(ray) || boundingBox.intersects(inverseRay))
		trees.push(this);
	else
		return;

	while(!trees.empty())
	{
		Octree * t = trees.top(); trees.pop();

		if (!t->intersectHit(tris))
		{
			for (StdVector<Octree>::iterator child = t->children.begin(); child != t->children.end(); child++)
			{
				if (child->boundingBox.intersects(ray) || child->boundingBox.intersects(inverseRay))
				{
					trees.push(&(*child));
				}
			}
		}
	}
}

//BaseTriangle* Octree::findClosestTri( const Ray & ray, IndexSet & tris, Mesh * mesh, HitResult & hitRes )
//{
//	double minDist = DBL_MAX;
//	BaseTriangle *closestFace = NULL, *curr_FaceType = NULL;

//	double u = 0.0, v = 0.0;
//	double actualMinDist = 0;

//	// Find the triangles in this tree
//	intersectRayBoth(ray, tris);

//	for(IndexSetIter it = tris.begin(); it != tris.end(); it++)
//	{
//		curr_FaceType = mesh->f(*it);

//		curr_FaceType->intersectionTest(ray, hitRes, true);

//		if(hitRes.hit && (abs(hitRes.distance) < minDist))
//		{
//			closestFace = curr_FaceType;

//			minDist = abs(hitRes.distance);
//			actualMinDist = hitRes.distance;

//			u = hitRes.u;
//			v = hitRes.v;
//		}
//	}

//	if(closestFace)
//	{
//		// set 'hitRes' to that closest hit result
//		hitRes.distance = actualMinDist;
//		hitRes.index = closestFace->index;
//		hitRes.u = u;	hitRes.v = v;
//		hitRes.hit = true;
//	}
//	else
//	{
//		hitRes.hit = false;
//		hitRes.distance = DBL_MAX;
//		hitRes.index = -1;
//	}

//	return closestFace;
//}

bool Octree::testIntersectHit( const Ray& ray, HitResult & hitRes )
{
	if(this->children.size() > 0)
		return false;

	// Do actual intersection test
	for(StdVector<BaseTriangle*>::iterator face = triangleData.begin(); face != triangleData.end(); face++)
	{
		BaseTriangle* f = *face;

		f->intersectionTest(ray, hitRes, true);

		if(hitRes.hit)
			return true;
	}

	return false;
}

void Octree::testIntersectRayBoth( const Ray& ray, HitResult & hitRes )
{
	stack<Octree*> trees;

	Ray inverseRay = ray.inverse();

	if (boundingBox.intersects(ray) || boundingBox.intersects(inverseRay))
		trees.push(this);

	while(!trees.empty())
	{
		Octree * t = trees.top(); trees.pop();

		if (t->testIntersectHit(ray, hitRes))
		{
			return;
		}
		else
		{
			for (StdVector<Octree>::iterator child = t->children.begin(); child != t->children.end(); child++)
			{
				if (child->boundingBox.intersects(ray) || child->boundingBox.intersects(inverseRay))
					trees.push(&(*child));
			}
		}
	}
}

StdVector<BaseTriangle*> Octree::getTriangleData()
{
	return triangleData;
}

Octree * Octree::root()
{
	if(parent == NULL)
		return this;
	else
		return parent->root();
}
