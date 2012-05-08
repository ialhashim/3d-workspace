#include "QuickMeshViewer.h"
#include <QFileInfo>

QuickMeshViewer::QuickMeshViewer( QWidget * parent /*= 0*/ ) :QGLViewer(parent)
{
	this->setMaximumSize(200,200);

	this->isActive = false;
	
	//connect(this, SIGNAL(meshLoaded()), SLOT(updateGL()));
}

void QuickMeshViewer::init()
{
	QGLViewer::init();

	setBackgroundColor(QColor(50,60,60));

	// Light
	GLfloat lightColor[] = {0.9f, 0.9f, 0.9f, 1.0f};
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor);

	// Camera
	camera()->setType(Camera::ORTHOGRAPHIC);

	resetView();
}

void QuickMeshViewer::draw()
{
	if(!this->isActive) return;

	glEnable(GL_MULTISAMPLE);
	mesh.draw();
}

void QuickMeshViewer::postDraw()
{
	if(!this->isActive) return;

	if(this->mesh.isLoading){
		glColor3d(1,1,1);
		renderText(8, 15, "Loading...");
	}
	else
	{
		glColor4d(1,1,1,0.5);
		renderText(8, 15, QFileInfo(mesh.fileName).baseName());
	}

	if(this->hasFocus())
	{
		int w = width(), h = height();
		startScreenCoordinatesSystem();
		glColor3d(0.8,0.2,0.2);
		glLineWidth(10);
		glBegin(GL_LINE_STRIP);
		glVertex3dv(Vec(0,0,0));
		glVertex3dv(Vec(w,0,0));
		glVertex3dv(Vec(w,h,0));
		glVertex3dv(Vec(0,h,0));
		glVertex3dv(Vec(0,0,0));
		glEnd();
		stopScreenCoordinatesSystem();
	}

	QGLViewer::postDraw();
}

void QuickMeshViewer::focusInEvent( QFocusEvent * event )
{
	emit(gotFocus(this));
}

void QuickMeshViewer::clearMesh()
{
	mesh.isLoading = true;
	mesh.clear();
}

void QuickMeshViewer::loadMesh( QString fileName )
{
	mesh.isLoading = true;
	mesh.load(fileName);
	mesh.isLoading = false;

	isActive = true;

	emit(meshLoaded());
}

void QuickMeshViewer::resetView()
{
	camera()->setSceneRadius(2.0);
	camera()->setUpVector(Vec(0,0,1));
	camera()->setSceneCenter(Vec(0,0,0));
	camera()->setPosition(Vec(1.25,1.25,1));
	camera()->lookAt(Vec(0,0,0));

	setGridIsDrawn(false);
	setAxisIsDrawn(false);
}
