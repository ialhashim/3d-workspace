#pragma once

#include <QObject>
#include "QSurfaceMesh.h"
#include "Voxeler.h"
#include "FFD.h"
#include "Graph.h"

class VoxelDeformer : public QObject{
	Q_OBJECT
public:
	VoxelDeformer(QSurfaceMesh * fromMesh, double voxel_size = 0.1);

	QSurfaceMesh * mesh;
	Voxeler * voxler;
	std::vector<Point> voxelPnts;				// control points positions
	std::vector< std::vector<double> > coord;	// mesh points coordinates
	
	std::vector<Point> meshPoints;
	std::vector<int> pointToFFD;

	StdVector<QControlPoint *> cpnts;

	StdVector<Point> startPos;

	double voxelSize;
	int voxelRange;

	std::vector<FFD*> ffd;

	Graph<int, double> voxelGraph;
	StdVector<double> computeNormalizedDistance(int startCornerIndex);

	QControlPoint * getQControlPoint( int selected );

	void draw();
	void drawNames();

	// debug
	std::vector<Vec3d> debugPoints, debugPoints2, debugPoints3;

	void select(int i);
	double sigmaControl;

private:
	int selectedCorner;
	
public slots:
	void update();
	void push(Vec3d from, Vec3d to);

signals:
	void meshDeformed();
};
