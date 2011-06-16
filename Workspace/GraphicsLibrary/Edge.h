#pragma once

#include "Point.h"
#include "Face.h"

class Edge
{
public:
	int vIndex[2];

	Edge(){ vIndex[0] = vIndex[1] = -1; }

	Edge(int v1, int v2)
	{
		vIndex[0] = v1;
		vIndex[1] = v2;
	}

	Edge(const Edge& newEdge)
	{
		vIndex[0] = newEdge.vIndex[0];
		vIndex[1] = newEdge.vIndex[1];
	}

	inline int operator[](int i)	const   { 	return vIndex[i];	}
	inline int neighbor()		const	{	return vIndex[1];	}

	float length(Vector<Point3D> & vertex)
	{
		return (vertex[vIndex[0]] - vertex[vIndex[1]]).norm();
	}
};

struct CompareEdge{
	bool operator()(const Edge &a, const Edge &b) 
	{
		if (a.vIndex[0] < b.vIndex[0]) return true;
		if (a.vIndex[0] > b.vIndex[0]) return false;
		return a.vIndex[1] < b.vIndex[1];
	}
};

typedef StdSet<Edge,CompareEdge> EdgeSet;
