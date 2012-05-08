#pragma once

#include "GraphicsLibrary/Mesh/QSurfaceMesh.h"
#include "RMF.h"

class GeneralizedCylinder{
private:
	QSurfaceMesh * src_mesh;

public:
	GeneralizedCylinder(){ src_mesh = NULL; }
	GeneralizedCylinder( std::vector<Point> spinePoints, QSurfaceMesh * mesh, bool computeRadius = true );
	QSurfaceMesh geometry;
	RMF frames;

	void realignCrossSections();

	void draw();
	bool isDrawFrames;

	// debug
	std::vector<Point> debugPoints;
	
public:
	class Circle { 
	public: 
		Point center; double radius; Normal n; uint index;

		Circle(Point newCenter, double newRadius, Normal newNormal, uint newIndex)
		{ center = newCenter; radius = newRadius; n = newNormal; index = newIndex; };

		std::vector<Point> toSegments( int numSegments, const Vec3d& startVec, double delta = 1.0);
		Normal normal();
	};

	std::vector<Circle> crossSection;

};
