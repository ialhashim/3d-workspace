#pragma once

// Qt stuff
#include <QColor>

// OpenMesh inclusion
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
using namespace OpenMesh;
typedef PolyMesh_ArrayKernelT<> OpenPolyMesh;

// Our interface
class QOpenMesh : public OpenPolyMesh
{
private:
	QColor color;

	// Bounding box
	Point bbox_min, bbox_max;

	void computeNormals();
	void computeBoundingBox();
	
public:
	void init();

	double radius();

	void draw();
};
