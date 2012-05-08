#pragma once 

#include "GraphicsLibrary/Mesh/QSurfaceMesh.h"

class Smoother{
private:
	static void treatBorders(Surface_mesh * m, std::vector<Surface_mesh::Vertex> & U);
	static void untreatBorders(Surface_mesh * m, std::vector<Surface_mesh::Vertex> & U);

public:
	static void LaplacianSmoothing(Surface_mesh * m, int numIteration, bool protectBorders = true);
	static Point LaplacianSmoothVertex(Surface_mesh * m, int vi);

	static void ScaleDependentSmoothing(Surface_mesh * m, int numIteration, double step_size = 0.5, bool protectBorders = true);
	static Point ScaleDependentSmoothVertex(Surface_mesh * m, Surface_mesh::Vertex v, double step_size = 0.5);

	static void MeanCurvatureFlow(QSurfaceMesh * m, int numIteration = 1, double step = 0.01, bool isVolumePreservation = true);
	static void MeanCurvatureFlowExplicit(QSurfaceMesh * m, int numIteration, double step = 0.5);
};
