#pragma once

#include "Scene.h"
#include <QMap>

class SubScene: public QObject{

	Q_OBJECT

public:
	SubScene( Scene * parentScene, int X, int Y, int Width, int Height );

	QString caption;

	Scene * activeScene;
	bool isSelected;
	bool isDraw;

	int x,y;
	int width,height;

	void drawFrame(Vec4d color = Vec4d(0,0,0,1));
	void draw();
	void postDraw();

	void select();

	QSegMesh * activeObject();

	bool contains(Vec2i pixel);

	// Auxiliary
	QMap<QString, double> property;
	QMap<QString, QVector<double>> properties;
	QMap<QString, int> property_int;
	QMap<QString, bool> property_bool;
};
