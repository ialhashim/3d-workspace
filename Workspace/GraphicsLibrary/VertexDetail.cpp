#include "VertexDetail.h"

VertexDetail::VertexDetail(int Index) : index(Index)
{
	ifaces = Vector<Face *> ();
	flag = VF_CLEAR;
}

StdSet<int> VertexDetail::adjacentVertices()
{
	StdSet<int> adjVertex;

	Face * face;

	for(Vector<Face *>::iterator f_it = ifaces.begin(); f_it != ifaces.end(); f_it++)
	{
		face = *f_it;

		for(int i = 0; i < 3; i++)
		{
			if(face->vIndex[i] != index)
				adjVertex.insert(face->vIndex[i]);
		}
	}

	return adjVertex;
}

bool VertexDetail::hasNeighbour(int vindex)
{
	StdSet<int> neighbours = adjacentVertices();

	return SET_HAS(neighbours, vindex);
}

bool VertexDetail::hasAnyNeighbour(StdSet<int> vindices)
{
	StdSet<int> neighbours = adjacentVertices();

	foreach(int i, neighbours)
	{
		foreach(int j, vindices)
		{
			if(i == j)
				return true;
		}
	}

	return false;
}

int VertexDetail::valence()
{
	return this->adjacentVertices().size();
}

bool VertexDetail::checkIsBorder()
{
        if(this->flag == VF_BORDER || this->valence() > (int)ifaces.size())
		return true;
	else
		return false;
}

bool VertexDetail::isMissingBorder()
{
        if(this->valence() > (int)ifaces.size() + 1)
		return true;
	else
		return false;
}

bool VertexDetail::checkMarkBorder()
{
        if(this->valence() > (int)ifaces.size())
	{
		this->flag = VF_BORDER;
		return true;
	}
	else
		return false;
}

bool VertexDetail::isBorderFlag()
{
	return this->flag == VF_BORDER;
}

void VertexDetail::insertFace(Face * newAdjFace)
{
	ifaces.push_back(newAdjFace);
}

void VertexDetail::detachFaces(const StdSet<Face *>& removeFaces)
{
	Vector<Face *> newAdj;

	foreach(Face * f, ifaces)
	{
		if(!SET_HAS(removeFaces, f))
		{
			newAdj.push_back(f);
		}
	}

	ifaces = newAdj;
}

void VertexDetail::attachFaces(const StdSet<Face *>& newFaces)
{
	foreach(Face * f, newFaces)
	{
		if(f->hasVertex(this->index))
			insertFace(f);
	}
}

void VertexDetail::unsetAndMark()
{
	ifaces.clear();
	flag = VF_DEAD;
}

Face * VertexDetail::sharedFace( int vi )
{
	Face * face = NULL;

	for(Vector<Face *>::iterator f_it = ifaces.begin(); f_it != ifaces.end(); f_it++)
	{
		face = *f_it;

		for(int i = 0; i < 3; i++)
		{
			if(face->vIndex[i] == vi)
				return face;
		}
	}

	return face;
}
