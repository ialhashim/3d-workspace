#include "global.h"
#include "Workspace.h"
#include <QVBoxLayout>
#include <QFileInfo>

Workspace::Workspace(QWidget *parent, Qt::WFlags flags)	: QMainWindow(parent, flags)
{
	ui.setupUi(this);

	QVBoxLayout * leftLayout = (QVBoxLayout *) ui.leftDockWidget->layout();
	QVBoxLayout * rightLayout = (QVBoxLayout *) ui.rightDockWidget->layout();

	tp = new TransformationPanel();
	rightLayout->addWidget(tp);
	
	mi = new MeshInfoPanel();
	rightLayout->addWidget(mi);

	// Connect to mesh management
	connect(ui.actionImportObject, SIGNAL(triggered()), mDoc, SLOT(importObject()));

	// Add new scene action
	connect(ui.actionNewScene, SIGNAL(triggered()), SLOT(addNewScene()));

	// Create new scene when we start by default
	sceneCount = 0;

	leftLayout->addStretch();
	rightLayout->addStretch();
}

Workspace::~Workspace()
{
	
}

void Workspace::addNewScene()
{
	Scene * newScene = new Scene(this);

	ui.sceneArea->addSubWindow(newScene);
	sceneCount++;

	newScene->showMaximized();
	newScene->setWindowTitle("Untitled");

	// Workspace window
	connect(newScene, SIGNAL(gotFocus(Scene*)), SLOT(setActiveScene(Scene*)));
	connect(newScene, SIGNAL(lostFocus(Scene*)), SLOT(disconnectScene(Scene*)));

	// MeshDoc
	connect(newScene, SIGNAL(objectDiscarded(QString)), mDoc, SLOT(deleteObject(QString)));

	// Object transformation
	connect(newScene, SIGNAL(gotFocus(Scene*)), tp, SLOT(setActiveScene(Scene*)));
	connect(tp, SIGNAL(objectModified()), newScene, SLOT(updateActiveObject()));

	connect(newScene, SIGNAL(gotFocus(Scene*)), mi, SLOT(setActiveScene(Scene*)));

	setActiveScene(newScene);
}

void Workspace::setActiveScene(Scene* scene)
{
	activeScene = scene;

	QString title = QString("%1 - %2")
		.arg(QFileInfo(QApplication::applicationFilePath()).baseName())
		.arg(scene->windowTitle());

	this->setWindowTitle(title);

	// View operations
	connect(ui.actionCameraProjection, SIGNAL(triggered()), activeScene, SLOT(toggleCameraProjection()), Qt::UniqueConnection);
	
	// Connect mDoc
	activeScene->connect(ui.actionExportObject, SIGNAL(triggered()), SLOT(exportActiveObject()), Qt::UniqueConnection);
	activeScene->connect(mDoc, SIGNAL(objectImported(QSegMesh*)), SLOT(setActiveObject(QSegMesh*)), Qt::UniqueConnection);
	activeScene->connect(mDoc, SIGNAL(printMessage(QString)), SLOT(print(QString)), Qt::UniqueConnection);
	mDoc->connect(activeScene, SIGNAL(exportActiveObject(QSegMesh*)), SLOT(exportObject(QSegMesh*)), Qt::UniqueConnection);

	// Set active scene
	tp->setActiveScene(activeScene);
	mi->setActiveScene(activeScene);
}

void Workspace::disconnectScene(Scene* scene)
{
	mDoc->disconnect();
	ui.actionCameraProjection->disconnect();
}
