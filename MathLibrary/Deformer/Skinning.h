#pragma once

#include "GraphicsLibrary/Mesh/QSurfaceMesh.h"
#include "GraphicsLibrary/Skeleton/GeneralizedCylinder.h"
#include <Eigen/Geometry>
using namespace Eigen;

class Skinning{

private:
	struct SkinningCoord{
		int		n1, n2;	// n1+1 = n2 : (n1, n2) is the closest segment
							// n1 == n2 : n1 is the closest node
		double	time;	// Time on the segment (n1, n2) or ( n1, n1+n1.norm() )
		Vec3d	d;		// Perpendicular offset from the skeleton

		SkinningCoord(){ }
		SkinningCoord(int i1, int i2, double t, Vec3d delta) : 
			n1(i1), n2(i2), time(t), d(delta) {}
	};

public:
	Skinning(QSurfaceMesh * src_mesh, GeneralizedCylinder * using_gc);

	void deform();
	std::vector<double> getCoordinate(Point p);
	Point fromCoordinates(std::vector<double> coords);

private:
	SkinningCoord	computeCoordinates(Point v);
	Point			fromCoordinates(SkinningCoord coords);
	void			computeMeshCoordinates();
	Matrix3d		rotationOfCurve(int cid);

private:
	QSurfaceMesh * mesh;
	GeneralizedCylinder * gc;
	GeneralizedCylinder origGC;
	std::vector< SkinningCoord > coordinates;
};
