#ifndef HALFEDGE_H
#define HALFEDGE_H

#include "Face.h"
#include "Edge.h"

class HalfEdge
{
public:
	Face* pFace[2];
	Edge edge;

	HalfEdge()
	{
		pFace[0] = pFace[1] = NULL;
	}

	HalfEdge(Edge & e, Face * f1, Face * f2)
	{
		edge.vIndex[0] = e.vIndex[0];
		edge.vIndex[1] = e.vIndex[1];

		pFace[0] = f1;
		pFace[1] = f2;
	}

	HalfEdge(const HalfEdge& other)
	{
		this->pFace[0] = other.pFace[0];
		this->pFace[1] = other.pFace[1];

		this->edge = other.edge;
	}

	inline Face * face(int i)	const { return pFace[i]; }
	inline int vertex(int i)	const { return edge.vIndex[i]; }
	
	Face * otherFace(Face * f) const 
	{ 
		if (pFace[0] == f) return pFace[1];
		if (pFace[1] == f) return pFace[0];

		return NULL;
	}

	bool isBorder() const
	{
		return !(pFace[0] && pFace[1]);
	}

	bool isNull() const
	{
		return (pFace[0] == NULL) && (pFace[1] == 0);
	}

	int fIndex(int i)
	{
		return pFace[i]->index;
	}

	Face * borderFace()
	{
		if(pFace[0] == NULL)
			return pFace[0];
		else
			return pFace[1];
	}

	inline bool hasVertex(int vi) const
	{
		return edge[0] == vi || edge[1] == vi;
	}

	bool hasEdge(const Edge & e) const
	{
		return hasVertex(e[0]) && hasVertex(e[1]);
	}

	bool hasFace(Face * f) const
	{
		return pFace[0] == f || pFace[1] == f;
	}
};

struct CompareHalfEdge{
	bool operator()(const HalfEdge &a, const HalfEdge &b) 
	{
		if (a.edge.vIndex[0] < b.edge.vIndex[0]) return true;
		if (a.edge.vIndex[0] > b.edge.vIndex[0]) return false;
		return a.edge.vIndex[1] < b.edge.vIndex[1];
	}
};

typedef StdSet<HalfEdge,CompareHalfEdge> HalfEdgeSet;
typedef Vector<HalfEdge> HalfEdgeVector;

#endif // HALFEDGE_H
