#pragma once

#include "GUI/Viewer/libQGLViewer/QGLViewer/qglviewer.h"
using namespace qglviewer;
#define GL_MULTISAMPLE 0x809D

#include "QuickMesh.h"

class LoaderThread;

class QuickMeshViewer : public QGLViewer{
	Q_OBJECT
public:
	QuickMeshViewer(QWidget * parent = 0);

	virtual void init();
	virtual void resetView();
	virtual void draw();
	virtual void postDraw();

	void focusInEvent( QFocusEvent * event );
	void clearMesh();

	bool isActive;
	QuickMesh mesh;

	QString meshFileName() { return mesh.fileName; }

public slots:
	void loadMesh(QString fileName);

signals:
	void gotFocus(QuickMeshViewer*);
	void meshLoaded();
};
