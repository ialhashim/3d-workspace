#pragma once

#include "Mesh.h"
#include "kdtree.h"
#include "Voxel.h"

#define glv glVertex3dv
#define gln glNormal3d

struct FaceBounds { int minX, minY, minZ; int maxX, maxY, maxZ;};

class Voxeler
{
private:
	Mesh * mesh;
	Vector<Voxel> voxels;
	KDTree kd;

	double voxelSize;

	Voxel minVox;
	Voxel maxVox;

	// Special voxels
	KDTree outerVoxels, innerVoxels;

public:
	Voxeler(Mesh * src_mesh = NULL, double voxel_size = 1.0);
	
	FaceBounds findFaceBounds( Face * f );
	bool isVoxelIntersects( const Voxel & v, Face * f );
	void computeBounds();

	// Find inside and outside of mesh surface
	Vector<Voxel> fillOther();
	void fillInsideOut(KDTree & inside, KDTree & outside);
	void fillOuter(KDTree & outside);

	// Visualization:
	void draw();
	void setupDraw();

	GLuint d1, d2;
};
