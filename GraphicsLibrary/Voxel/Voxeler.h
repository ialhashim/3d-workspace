#pragma once

#include "GraphicsLibrary/Mesh/QSurfaceMesh.h"
#include "GraphicsLibrary/SpacePartition/kdtree.h"
#include "Voxel.h"
#include "MathLibrary/Bounding/BoundingBox.h"

#define glv glVertex3dv
#define gln glNormal3d

class Voxeler
{
private:
	QSurfaceMesh * mesh;
	KDTree kd;

	// Special voxels
	KDTree outerVoxels, innerVoxels;

public:
	Voxeler( QSurfaceMesh * src_mesh, double voxel_size, bool verbose = false);

	FaceBounds findFaceBounds( Surface_mesh::Face f );
	bool isVoxelIntersects( const Voxel & v, Surface_mesh::Face f );
	
	void update();
	void computeBounds();

	// Grow larger by one voxel
	void grow();

	// Find inside and outside of mesh surface
	std::vector< Voxel > fillOther();
	void fillInsideOut(KDTree & inside, KDTree & outside);
	void fillOuter(KDTree & outside);

	// Intersection
	std::vector<Voxel> Intersects(Voxeler * other);
	std::map<int, Voxel> around(Point p);

	// Visualization:
	void draw();
	void setupDraw();
	static void drawVoxels( const std::vector< Voxel > & voxels, double voxel_size = 1.0);

	KDTree corner_kd;
	std::vector< Point > corners;
	std::vector< std::vector<int> > cornerIndices;
	std::vector< int > cornerCorrespond;
	std::vector< Point > getCorners(int vid);
	int getClosestVoxel(Vec3d point);
	int getEnclosingVoxel( Vec3d point );

	std::vector< Point > getVoxelCenters();

	std::vector< Voxel > voxels;
	int getVoxelIndex(Voxel v);
	uint d1, d2;
	bool isVerbose;
	bool isReadyDraw;

	double voxelSize;

	Voxel minVox;
	Voxel maxVox;
};
