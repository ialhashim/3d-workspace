#pragma once

#include "ui_RotationWidget.h"
#include "GUI/Scene.h"

class TransformationPanel : public QWidget
{
	Q_OBJECT

public:
	TransformationPanel();

	Ui::RotationWidget rotWidget;

	Scene * activeScene;
	QSegMesh* activeObject();

public slots:
	void setActiveScene( Scene * scene );
	void applyRotation();

signals:
	void objectModified();
};
