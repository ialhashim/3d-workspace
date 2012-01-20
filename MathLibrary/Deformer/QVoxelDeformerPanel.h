#pragma once

#include "ui_QVoxelDeformerWidget.h"
#include "Scene.h"
#include "VoxelDeformer.h"

class QVoxelDeformerPanel : public QWidget
{
	Q_OBJECT

private:
	Ui::QVoxelDeformerWidget dw;
	Scene * activeScene;

public:
	QVoxelDeformerPanel();

	VoxelDeformer * activeDeformer;

public slots:
	void setActiveScene(Scene *);
	void onBuildButtonClicked();
	void onFalloffChanged(double);

signals:
	void deformerCreated( VoxelDeformer * );
	void updateScene();

};
