#include "TransformationPanel.h"

TransformationPanel::TransformationPanel()
{
	rotWidget.setupUi(this);
	
	connect(rotWidget.rotateButton, SIGNAL(clicked()), SLOT(applyRotation()));

	activeScene = NULL;
}

void TransformationPanel::setActiveScene( Scene * scene )
{
	if(activeScene != scene)
		activeScene = scene;
}

QSegMesh* TransformationPanel::activeObject()
{
	if (activeScene)
		return activeScene->activeObject();
	else 
		return NULL;
}

void TransformationPanel::applyRotation()
{
	QSegMesh * mesh = activeObject();

	double xAngle = RADIANS(rotWidget.xAngle->value());
	double yAngle = RADIANS(rotWidget.yAngle->value());
	double zAngle = RADIANS(rotWidget.zAngle->value());

	foreach(QSurfaceMesh * seg, mesh->getSegments())
	{
		qglviewer::Quaternion qx(qglviewer::Vec(1,0,0), xAngle);
		qglviewer::Quaternion qy(qglviewer::Vec(0,1,0), yAngle);
		qglviewer::Quaternion qz(qglviewer::Vec(0,0,1), zAngle);

		Surface_mesh::Vertex_property<Point> points = seg->vertex_property<Point>("v:point");
		Surface_mesh::Vertex_iterator vit, vend = seg->vertices_end();

		for (vit = seg->vertices_begin(); vit != vend; ++vit){
			qglviewer::Vec v(points[vit].x(), points[vit].y(), points[vit].z());
			v = (qx * qy * qz).rotate(v);
			points[vit] = Point(v.x, v.y, v.z);
		}
	}

	// Reset to zero
	//rotWidget.xAngle->setValue(0);
	//rotWidget.yAngle->setValue(0);
	//rotWidget.zAngle->setValue(0);

	// Tell blah to update
	emit(objectModified());
}
