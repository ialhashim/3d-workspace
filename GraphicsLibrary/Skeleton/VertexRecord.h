#pragma once

#include <float.h>
#include "GraphicsLibrary/Mesh/QSurfaceMesh.h"

// data structures
#include <queue>
#include <QList>
#include <QMap>
#include <QSet>
#include <QVector>

#include <Eigen/Geometry>
using namespace Eigen;

class VertexRecord {

public:
	Vec3d pos;

	std::set<uint> adjV;
	std::set<uint> adjF;
	std::set<uint> collapseFrom;

	std::map<uint, double> boundaryLength;

	int vIndex;
	int minIndex;
	double minError;
	bool center;
	double err;
	double nodeSize;
	int PQIndex;

	MatrixXd matrix;

	VertexRecord()
	{
		// Defaults
		vIndex = minIndex = -1;
		minError = std::numeric_limits<double>::max();
		center = false;
		err = 0;
		nodeSize = 1;
		matrix = Matrix4d::Zero(4,4);
		PQIndex = -1;
	}

	VertexRecord(QSurfaceMesh & mesh, const Surface_mesh::Vertex & v)
	{
		// Defaults
		vIndex = minIndex = -1;
		minError = std::numeric_limits<double>::max();
		center = false;
		err = 0;
		nodeSize = 1;
		matrix = Matrix4d::Zero(4,4);
		PQIndex = -99;

		Surface_mesh::Vertex_property< std::set<uint> > adjVV = mesh.vertex_property< std::set<uint> >("v:adjVV");
		Surface_mesh::Vertex_property< std::set<uint> > adjVF = mesh.vertex_property< std::set<uint> >("v:adjVF");

		vIndex = v.idx();

		pos = mesh.getVertexPos(v);

		adjV = adjVV[v];
		adjF = adjVF[v];
	}

	int CompareTo(const VertexRecord &other) const
	{ 
		if (minError < other.minError) return -1;
		if (minError > other.minError) return 1;
		return 0;
	}
};
