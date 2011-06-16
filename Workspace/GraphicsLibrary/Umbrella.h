#ifndef UMBRELLA_H
#define UMBRELLA_H

#include "VertexDetail.h"
#include "HalfEdge.h"

class Umbrella
{
public:
	int index;
	int flag;

	Vector<Face *> ifaces;
	HalfEdgeSet halfEdge;
	Vector<int> neighbor;

	int valence();

	Umbrella(VertexDetail * vd);
	Umbrella(){ index = -1; };

	Umbrella(const Umbrella& from)
	{  
		this->index = from.index;
		this->flag = from.flag;
		this->ifaces = from.ifaces;
		this->halfEdge = from.halfEdge;
		this->neighbor = from.neighbor;
	}

	Umbrella& operator= (const Umbrella& from)
	{  
		if (this != &from) {
			this->index = from.index;
			this->flag = from.flag;
			this->ifaces = from.ifaces;
			this->halfEdge = from.halfEdge;
			this->neighbor = from.neighbor;
		}
		return *this;
	}

	// Find edges and neighbors
	void loadHalfEdgeSet();

	// Border operations
	PairInt borderNeighbours();
	PairFaces borderFaces();

	PairFaces sharedFaces(Umbrella * other);
	HalfEdge sharedEdge(Umbrella * other);

	bool hasNeighbour(int index);
	bool hasNeighbour(Umbrella * other);

	bool shareTriangle(Umbrella * b, Umbrella * c);

	StdSet<int> faceNeighbours(Face * face);
	PairFaces faceNeighboursPair(Face * face);

	inline Vec vec() { return *ifaces.back()->PointIndexed(index); }

	Normal normal();

	// Scale-dependent umbrella operator
	Vec scaleDependentFaired();
};

#endif // UMBRELLA_H

