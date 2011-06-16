#ifndef SLICER_H
#define SLICER_H

#include "Mesh.h"
#include "Plane.h"

#include "HashTable.h"

struct SliceResult
{
	StdSet<Face*> cutFaces, removeFaces, keepFaces;
	IntSet cutPoints;
	IntSet newPoints;
	EdgeSet cutEdges;

	SliceResult(StdSet<Face*>& faces, IntSet& points, IntSet& new_points, 
		EdgeSet& edges, StdSet<Face*>& RemoveFaces, StdSet<Face*>& KeepFaces)
	{
		cutPoints = points; 
		newPoints = new_points;

		cutFaces = faces; 
		cutEdges = edges; 

		removeFaces = RemoveFaces; 
		keepFaces = KeepFaces;
	}

	SliceResult(){}
};

class Slicer
{
private:

public:
	static void FindSliceStrip(Mesh * mesh, Vector<int> & facesIndices, const Plane& cutPlane, 
							StdSet<int> & cutPoints, IntSet & cutPointsTable, IntSet & cutFaces, EdgeSet & cutEdges);
        static SliceResult SliceAt(Mesh * mesh, Vector<int> & facesIndices,  const Plane& cutPlane);
};

#endif // SLICER_H
