#pragma once

#include <QColorDialog>
#include <QQueue>

#include "QGLViewer/qglviewer.h"
using namespace qglviewer;

#include "QMesh.h"
#include "QOpenMesh.h"

enum ViewMode { VIEW, SELECTION, MODIFY };
enum SelectMode { NONE, MESH, SKELETON_NODE, SKELETON_EDGE, SKELETON_FACES, RECONSTRUCTED_POINTS, VERTEX};
enum ModifyMode { DEFAULT, CP_REF_VECTOR, MOVE_VERTEX };

class Scene : public QGLViewer{

	Q_OBJECT

public:
	Scene(QWidget *parent = 0);

	// Setup scene
	virtual void init();
	void setupCamera();
	void setupLights();

	// OpenGL Drawing
	virtual void draw();
	virtual void drawWithNames();
	virtual void postDraw();

	// Scene Visualizations
	void drawCornerAxis();

	// Mouse & Keyboard stuff
	virtual void mousePressEvent(QMouseEvent* e);
	virtual void mouseReleaseEvent(QMouseEvent* e);
	virtual void mouseMoveEvent(QMouseEvent* e);
	virtual void keyPressEvent(QKeyEvent *e);

	// SELECTION
	virtual void postSelection(const QPoint& point);

	// STATE
	ViewMode viewMode;
	SelectMode selectMode;
	ModifyMode modifyMode;

	QColor backColor;

	void setViewMode(ViewMode toMode);
	void setSelectMode(SelectMode toMode);
	void setModifyMode(ModifyMode toMode);

// TEXT ON SCREEN
public slots:
	void print(QString message, long age = 1000);
	void dequeueLastMessage();

private:
	QQueue<QString> osdMessages;
	QTimer *timer;

// Objects in the scene
private:
	QMap<QString, QMesh *> objects;
	QMap<QString, QOpenMesh *> openObjects;

public slots:
	void insertObject( QString fileName );
};
