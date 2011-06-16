#ifndef VERTEXDETAIL_H
#define VERTEXDETAIL_H

#include "Face.h"

enum VertexFlag{
    VF_CLEAR		=   1,
    VF_VISISTED		=   2,
    VF_SELECTED		=   4,
    VF_BORDER		=   8,
    VF_DEAD			=  16,
	VF_MODIFIED		=  32,
	VF_BAD_BORDER   =  64,
};

class VertexDetail
{
public:
	int index;
	Vector<Face *> ifaces;
	int flag;

	VertexDetail(int Index = -1);

	int valence();

	bool checkIsBorder();
	bool checkMarkBorder();
	bool isBorderFlag(); // simply check flag

	bool isMissingBorder();

	StdSet<int> adjacentVertices();

	bool hasNeighbour(int vindex);
	bool hasAnyNeighbour(StdSet<int> vindices);

	void insertFace(Face * newAdjFace);

	void detachFaces(const StdSet<Face *>& removeFaces);
	void attachFaces(const StdSet<Face *>& newFaces);

	void unsetAndMark();

	Face * sharedFace(int vi);
};


#endif // VERTEXDETAIL_H
