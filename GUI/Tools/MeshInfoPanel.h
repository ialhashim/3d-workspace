#pragma once

#include <QWidget>
#include "GUI/Scene.h"

#include "ui_MeshInfo.h"

class MeshInfoPanel : public QWidget
{
	Q_OBJECT

public:
	MeshInfoPanel();

	Ui::MeshInfoWidget infoWidget;

	Scene * activeScene;
	QSegMesh* activeObject();

public slots:
	void setActiveScene( Scene * scene );

signals:
	void objectModified();
};
