#pragma once

#include "ui_DeformerWidget.h"
#include "GUI/Scene.h"
#include "QFFD.h"

class DeformerPanel : public QWidget
{
	Q_OBJECT

private:
	Ui::DeformerWidget dw;
	Scene * activeScene;

public:
	DeformerPanel();

	QFFD * activeDeformer;

public slots:
	void setActiveScene(Scene *);
	void onCreateBoundingClicked();

signals:
	void deformerCreated( QFFD * );
};
