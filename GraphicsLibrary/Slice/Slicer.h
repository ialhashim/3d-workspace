#pragma once

#include "QSurfaceMesh.h"
#include "Plane.h"

#include "HashTable.h"

struct SliceResult
{
	StdSet<Face*> cutFaces, removeFaces, keepFaces;
	IntSet cutPoints;
	IntSet newPoints;
	EdgeSet cutEdges;

	SliceResult(StdSet<Surface_mesh::Face>& faces, IntSet& points, IntSet& new_points, 
		EdgeSet& edges, StdSet<Surface_mesh::Face>& RemoveFaces, StdSet<Surface_mesh::Face>& KeepFaces)
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
	static void FindSliceStrip(QSurfaceMesh * mesh, StdVector<int> & facesIndices, const Plane& cutPlane, 
							StdSet<int> & cutPoints, IntSet & cutPointsTable, IntSet & cutFaces, EdgeSet & cutEdges);
        static SliceResult SliceAt(QSurfaceMesh * mesh, StdVector<int> & facesIndices,  const Plane& cutPlane);
};
