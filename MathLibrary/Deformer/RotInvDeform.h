#pragma once

#include "GraphicsLibrary/Mesh/QSurfaceMesh.h"
#include "Frame3.h"

typedef Vec3d Vector3;

class RotInvDeform{

public:
	RotInvDeform(QSurfaceMesh * usingMesh);

	QSurfaceMesh * mesh;

private:
	Surface_mesh::Vertex_property<Point> points;
	Surface_mesh::Vertex_property<Normal> normals;
	Surface_mesh::Vertex_iterator vit, vend;
	int nVerts;

	// Vertex info
	struct VtxInfo {
		Surface_mesh::Vertex vID;
		std::vector<Surface_mesh::Vertex> vNbrs;	// Neighbors 
		std::vector<double> vNbrWeights;			// Neighbors weights
		double fVertexArea;
		double fVertexWeight;
		double fVertexScale;

		Surface_mesh::Vertex nTangentNbr;	// nbr edge that is used to compute x axis of tangent frame
		Frame3 vFrame;						// original tangent frame
		Vector3 vFrameLaplacian;			// laplacian vector encoded in original frame

		Frame3 vTransFrame;					// current frame (solved for)
		Vector3 vLaplacian;					// current laplacian

		VtxInfo() { }
	};
	std::vector<VtxInfo> vVertices;

	// Rotational constraint
	struct RotConstraint {
		Surface_mesh::Vertex vID;
		Frame3 vFrame;
		double fWeight;
	};
	std::vector<RotConstraint> vRotConstraints;

	// Positional constraint
	struct PosConstraint {
		Surface_mesh::Vertex vID;
		Vector3 vPosition;
		double fWeight;
	};
	std::vector<PosConstraint> vPosConstraints;

public:

	void AddBoundaryConstraints(double fWeight = 1.0);

	void UpdatePositionConstraint( Surface_mesh::Vertex vID, const Vector3 & vPosition, double fWeight );
	void UpdateOrientationConstraint( Surface_mesh::Vertex vID, const Frame3 & vFrame, double fWeight );

private:

	Frame3 GetCurrentFrame( Surface_mesh::Vertex vID );
	Frame3 GetSolvedFrame( Surface_mesh::Vertex vID );

	void SetVertexWeight( Surface_mesh::Vertex vID, double fWeight );
	void SetVertexScale( Surface_mesh::Vertex vID, double fScale );

private:

	int nEdges;
	double avgVtxArea;
	double globalScale;
	bool isMatrixValid;
	
	void ComputeWeights();

	// Utility
	void CotangentWeights( Surface_mesh::Vertex v, std::vector<double> & vWeights, bool bnormalize = false );
	double VectorCot( const Vector3 & v1, const Vector3 & v2 );
	Vector3 MeshLaplacian( Surface_mesh::Vertex v, std::vector<double> & vWeights );

	// Setup systems
	void UpdateMatrices();

	// Orientation
	void UpdateRotMatrix();
	SparseMatrix<double,RowMajor> RotA;		
	MatrixXd RHSRot;
	void UpdateRHSRot();

	// Position
	void UpdatePosMatrix();
	SparseMatrix<double,RowMajor> PosA;
	MatrixXd RHSPos;
	void UpdateRHSPos();

public:
	void Solve();

	double GetLaplacianError();

	void draw();
};
