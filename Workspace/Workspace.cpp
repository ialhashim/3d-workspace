#include "Workspace.h"

#include "Scene.h"

Workspace::Workspace(QWidget *parent, Qt::WFlags flags)
	: QMainWindow(parent, flags)
{
	ui.setupUi(this);

	// New scene action
	connect(ui.actionNewScene, SIGNAL(triggered()), SLOT(addNewScene()));
	connect(ui.actionImportObject, SIGNAL(triggered()), SLOT(importObject()));

	// Create new scene when we start by default
	addNewScene();
}

Workspace::~Workspace()
{

}

void Workspace::addNewScene()
{
	Scene * newScene = new Scene;

	ui.sceneArea->addSubWindow(newScene);

	newScene->showMaximized();
}

void Workspace::importObject()
{
	Scene * selectedScene = static_cast<Scene*>(ui.sceneArea->activeSubWindow()->widget());

	QString fileName = QFileDialog::getOpenFileName(this, tr("Insert Mesh"), "", tr("Mesh Files (*.obj *.off)"));

	if(fileName.length())
		selectedScene->insertObject(fileName);

	emit(importedObject(fileName));
}
