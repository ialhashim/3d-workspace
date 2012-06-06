// Adapted from code by Shaoting Zhang
// http://sourceforge.net/projects/meshtools/
#pragma once

#include "GraphicsLibrary/Mesh/SurfaceMesh/Surface_mesh.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <Eigen/Geometry>
using namespace Eigen;

class ARAPDeformer{

public:
	ARAPDeformer(Surface_mesh * usingMesh);
	~ARAPDeformer();

private:
	void ComputeCotWeights();
	void BuildAndFactor();
	void SVDRotation();

public:
	void Deform(int ARAPIteration = 1);

private:
	std::vector<Matrix3d> R;
	std::vector<Vector3d> OrigMesh;
	std::vector<VectorXd> xyz;
	std::vector<VectorXd> b;

	// Frequently used
	int nVerts;
	Surface_mesh::Vertex_property<Point> points;
	Surface_mesh::Vertex_property<Normal> normals;
	Surface_mesh::Vertex_iterator vit, vend;
	Surface_mesh::Vertex_property< std::map<Surface_mesh::Vertex, double> > wij_weight;
	Surface_mesh::Vertex_property< bool > isAnchorPoint, isControlPoint;
	Surface_mesh::Vertex_around_vertex_circulator vvit, vvend;

	SparseMatrix<double> At;
	SimplicialLLT< SparseMatrix<double> > solver;
	bool isSolverReady;

public:
	Surface_mesh * mesh;

	// Control points
	void SetAnchor( const Surface_mesh::Vertex & v ){ isAnchorPoint[v] = true; }

	void UpdateControl( const Surface_mesh::Vertex & v, const Vec3d & newPos ){
		isControlPoint[v] = true;
		points[v] = newPos;
		isSolverReady = false;
	}

	void ClearAnchors(){
		for (vit = mesh->vertices_begin(); vit != vend; ++vit)
			isAnchorPoint[vit] = false;
		isSolverReady = false;
	}

	void ClearControl(){
		for (vit = mesh->vertices_begin(); vit != vend; ++vit)
			isControlPoint[vit] = false;
		isSolverReady = false;
	}

	void ClearAll(){
		ClearAnchors();
		ClearControl();
	}
};
