#pragma once

#include "GraphicsLibrary/Mesh/QSurfaceMesh.h"
#include "GUI/Viewer/libQGLViewer/QGLViewer/manipulatedFrame.h"

// Movable control points
class QControlPoint : public qglviewer::ManipulatedFrame
{
	Q_OBJECT

public:
	Vec3d pos;
	double weight;
	uint idx;
	Vec3i gridIdx;

	QControlPoint(const Vec3d & fromPos, uint index = 0, const Vec3i & gridIndex = Vec3i(0,0,0), double fromWeight = 1.0){
		pos = fromPos;
		idx = index;
		gridIdx = gridIndex;
		weight = fromWeight;

		setPosition(qglviewer::Vec(pos.x(), pos.y(), pos.z()));
		setConstraint(NULL);

		connect(this, SIGNAL(manipulated()), SLOT(updatePosition()));
	}
	
public slots:
	void updatePosition() {	
		pos.x() = this->position().x;
		pos.y() = this->position().y;
		pos.z() = this->position().z;
	}
};
