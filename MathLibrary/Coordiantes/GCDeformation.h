#pragma once

#include "GraphicsLibrary/Mesh/QSurfaceMesh.h"

class GCDeformation{

public:
	GCDeformation(QSurfaceMesh * forShape, QSurfaceMesh * usingCage);
	
	void deform();

	QSurfaceMesh * shape;
	QSurfaceMesh * cage;

	std::vector<Point> orginalCagePos, deformedCagePos;
	std::vector<Normal> orginalCageNormal, deformedCageNormal;
	std::vector<double> S;

	struct GreenCoordiante{
		std::vector<double> coord_v;
		std::vector<double> coord_n;
		bool insideCage;
		bool valid;
	};

	std::vector< GreenCoordiante > coords;

	void initDeform();

	Point deformedPoint(GreenCoordiante gc);

	GCDeformation::GreenCoordiante computeCoordinates(Vec3d point);
private:
	double GCTriInt(const Vec3d& p, const Vec3d& v1, const Vec3d& v2, const Vec3d& e);
};
