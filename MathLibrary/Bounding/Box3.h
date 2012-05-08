#pragma once

#include "GraphicsLibrary/Mesh/QSurfaceMesh.h"

class Box3
{
public:
	// Constructor
	Box3();
	Box3( Box3 &box );

	// Regularize
	void normalizeAxis();
	void makeRightHanded();
	void sort();

	// Operator
	bool operator == (Box3& box);
 

	// Proximity
	Vec3d ClosestPoint(const Vec3d& p);
	void ClosestSegment( Box3 other, Vec3d & p, Vec3d & q);
	Vec3d ClosestAxis( const Vec3d& v );

public:
	Point Center;
	std::vector<Vec3d> Axis;
	Vec3d Extent;

	std::vector<double> faceScaling;
};
