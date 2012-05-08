#pragma once

#include "GraphicsLibrary/Mesh/QSurfaceMesh.h"
#include "QControlPoint.h"

enum FFD_FitType {BoundingBoxFFD, VolumeFFD};

class FFD
{
public:
	FFD(QSurfaceMesh * src_mesh = NULL, FFD_FitType fit_type = BoundingBoxFFD);
	QSurfaceMesh * mesh;
	double width, length, height;

	Vec3i resolution;

	StdVector<QControlPoint*> points;
	StdVector< StdVector< StdVector < int > > > pointsGridIdx; 

	void bbFit(Vec3i res);

	void apply();

	Vec3d deformVertexLocal( const Vec3d & localPoint );
	Vec3d getWorldCoordinate(const Vec3d & pLocal);
	Vec3d getLocalCoordinates( const Vec3d & p );
	Vec3d mP, mS, mT, mU;	// the local frame coordinates
	StdVector<Vec3d> meshVerticesLocal;

	// Manually setup FFD
	void fixed( Vec3i res, Vec3d location, double spacing, StdMap<int,Point> pnts );
	StdMap<int,Point> applyFixed();
	StdMap<int,Point> fixedPointsLocal;

	// Debug
	StdVector<Vec3d> dbPoints;
	StdVector< Pair<Vec3d,Vec3d> > dbLines;
};

// math helpers
inline int fact(int x){
	if ( x == 0 ) return 1;
	int f = 1;
	for (int i=1;i<=x;++i) f*=i;
	return f;
}
inline int C(int a, int b){
	// return the aCb (binomial coefficients)
	assert(a>=b);
	return fact(a) / (fact(b) * fact(a-b));
}
