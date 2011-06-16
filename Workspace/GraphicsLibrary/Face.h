#pragma once

#include "Point.h"
#include "Edge.h"
#include "Triangle.h"

enum FaceFlag{
	FF_CLEAR =  -1,
	FF_INVALID_VINDEX = 2
};

struct VertexAngle{
	int index;
	double angle;
	VertexAngle(int vertIndex, double angleAt) : index(vertIndex), angle(angleAt){};
};

class Face : public BaseTriangle
{
public:
	Point3D * v[3];
	int vIndex[3];

	Face(){v[0] = v[1] = v[2] = NULL; vIndex[0] = vIndex[1] = vIndex[2] = -1; index = -1; flag = -1;}
	Face(int vIndex1, int vIndex2, int vIndex3, Point3D * v1, Point3D * v2, Point3D * v3, int Index);
	Face(const Face& from);
	Face& operator= (const Face& from);

	// DIRECT VERTEX ACCESS
	inline Point3D * operator[]( int subscript ) const { return v[ subscript ]; }
	inline Point3D * & operator[]( int subscript ) { return v[ subscript ]; }
	inline Point3D * P(int subscript) { return v [ subscript ]; }

	// VERTEX ACCESS
	inline int VIndex(const int & i){ return vIndex[i]; }
	Point3D * PointIndexed(int i);
	int closestVertexIndex(Point3D * vertex);
	int LocalIndexPointIndexed(int vi);
	int oppositeVertex(Edge & e);
	int otherVertexIndex(int vi);
	int otherVertexIndex(int vi1, int vi2);
	bool hasVertex(int i);
	void replacePoint(int oldIndex, int newIndex);
	void unsetVertexByIndex(int vi);
	void unset();

	Vec vec(int i) const;

	// COMPUTE NORMAL
	Vec normal() const;

	// COMPUTE CENTER & POINT from Barycentric
	Vec center() const;
	Vec getBary(double U, double V) const;

	// EDGE OPERATIONS
	bool hasEdge(const Edge & e);
	bool sharesEdge(Face & other);
	int sharesEdge(Face * other, int firstVertex);
	const EdgeSet edgesWithVertex(const int & vertex);
	double longestEdgeLength();
	Edge oppositeEdge(int vi);
	bool inOrder(const Edge & e);

	// AREA & VOLUME
	double area();
	double volume();

	// ANGLES
	double angle(int a);
	double angleCotangent(int a);
	double angleCotangent(Edge & e);
	VertexAngle largestAngle();

	Pair<Vec, Vec> spanAt(int a);

	// INTERSECTION FUNCTIONS
	void intersectionTest(const Ray & ray, HitResult & res, bool allowBack = false) const;
	void intersectionTest2(const Ray & ray, HitResult & res);

	double pointToPointsDistance(Vec & p);

	StdString toString();
};

#define PairFaces Pair<Face *, Face *>
