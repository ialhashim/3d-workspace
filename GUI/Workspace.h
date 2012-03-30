#pragma once

#include <QtGui/QMainWindow>
#include "ui_Workspace.h"

#include <QMdiSubWindow>

#include "QMeshDoc.h"
#include "Scene.h"
#include "GUI/Tools/TransformationPanel.h"
#include "GUI/Tools/MeshInfoPanel.h"

class Workspace : public QMainWindow
{
	Q_OBJECT

public:
	Workspace(QWidget *parent = 0, Qt::WFlags flags = 0);
	~Workspace();

	Ui::WorkspaceClass ui;
	Scene * activeScene;

public slots:
	void addNewScene();
	void setActiveScene(Scene* scene);
	void disconnectScene(Scene* scene);

private:
	MeshInfoPanel * mi;
	TransformationPanel * tp;
	
	int sceneCount;
};
