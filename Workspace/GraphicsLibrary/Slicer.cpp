#include "Slicer.h"

void Slicer::FindSliceStrip(Mesh * mesh, Vector<int> & facesIndices, 
							const Plane& plane, 
							StdSet<int> & cutPoints, IntSet & cutPointsTable, 
							IntSet & cutFaces, EdgeSet & cutEdges)
{
	// POINTS
	StdSet<int> cutPointsFirst;

	// FACES
	IntSet cutFacesFirst;

	// EDGES
	EdgeSet cutEdgesFirst;

	// Half Edges
	Vector<Umbrella> U;
	mesh->getUmbrellas(U);

	// First pass test
        for(int i = 0; i < (int)facesIndices.size(); i++)
	{
		Face * face = mesh->f(facesIndices[i]);

		int v0Index = face->vIndex[0], v1Index = face->vIndex[1], v2Index = face->vIndex[2];

		Vec v0 = face->vec(0);
		Vec v1 = face->vec(1);
		Vec v2 = face->vec(2);

		int vtest0 = 0, vtest1 = 0, vtest2 = 0;

		if(plane.IsFront(v0)) vtest0 = 1;
		if(plane.IsFront(v1)) vtest1 = 1;
		if(plane.IsFront(v2)) vtest2 = 1;

		int test = vtest0 + vtest1 + vtest2;

		if(test == 1 || test == 2)
		{
			if(vtest0) cutPointsFirst.insert(v0Index);
			if(vtest1) cutPointsFirst.insert(v1Index);
			if(vtest2) cutPointsFirst.insert(v2Index);

			if(vtest0 && vtest1) cutEdgesFirst.insert(Edge(v0Index, v1Index));
			if(vtest1 && vtest2) cutEdgesFirst.insert(Edge(v1Index, v2Index));
			if(vtest2 && vtest0) cutEdgesFirst.insert(Edge(v2Index, v0Index));

			cutFacesFirst.insert(face->index);
		}
	}

	// Second pass point tests
	foreach(const int& i, cutPointsFirst)
	{
		// Check if its surounded by cut faces, if so it is a redundant cut point
		int numFaces = 0;
		foreach(Face * face, U[i].ifaces)
		{
			if(cutFacesFirst.has(face->index))
				numFaces++;
		}

		// Add points that are useful
                if(numFaces < (int)U[i].ifaces.size())
		{
			cutPoints.insert(i);
			cutPointsTable.insert(i);
		}
	}

	// Make sure its a strip
	foreach(const int& i, cutPoints)
	{
		foreach(Face * face, U[i].ifaces)
		{
			if(cutFacesFirst.has(face->index))
			{
				cutFaces.insert(face->index);
			}
		}
	}

	// Filter out useless edges
	foreach(Edge e, cutEdgesFirst)
	{
		if(cutPointsTable.has(e[0]) && cutPointsTable.has(e[1]))
		{
			cutEdges.insert(e);

			//Debug edges:
			//mesh->testLines.push_back(Line(*mesh->v(e[0]), *mesh->v(e[1])));
		}
	}
}

SliceResult Slicer::SliceAt(Mesh * mesh, Vector<int> & facesIndices,  const Plane& cutPlane)
{
	StdSet<int> cutPoints;
	IntSet cutPointsTable, cutFaces, newPoints;
	EdgeSet cutEdges;

	// debug: draw slice plane
	//mesh->testPlanes.push_back(cutPlane);

	FindSliceStrip(mesh, facesIndices, cutPlane, cutPoints, cutPointsTable, cutFaces, cutEdges);

	// New vertices
	int vIndex = mesh->vertex.size();

	StdSet<Face *> modifiedFaces, allRemoveFaces, allKeepFaces;

	//DEBUG faces
	//mesh->greenFaces.push_back(StdSet<Face*>());
	//mesh->yellowFaces.push_back(StdSet<Face*>());

	foreach(const int& i, cutPoints)
	{
		// Creating new vertex
		mesh->addVertex(*mesh->v(i), vIndex++);
		mesh->vColor.push_back(mesh->vColor[i]);
		mesh->vNormal.push_back(Vec());

		// Get details for old and new
		VertexDetail * prev_vd = mesh->vd(i);
		VertexDetail * new_vd = &mesh->vertexInfo.back();

		newPoints.insert(new_vd->index);

		// debug: change color of old
		//mesh->vColor[prev_vd->index] = Color4(0,120,250);
		//mesh->testVertex.push_back(*mesh->v(prev_vd->index));

		Vector<Face *> prev_ifaces = prev_vd->ifaces;
		StdSet<Face *> removeFaces, keepFaces;

		// Find faces need to be removed from old vertex
		foreach(Face * f, prev_ifaces)
		{
			if(cutFaces.has(f->index))
			{
				f->replacePoint(prev_vd->index, new_vd->index);

				removeFaces.insert(f);
				allRemoveFaces.insert(f);

				// DEBUG:
				//mesh->yellowFaces.back().insert(f);
			}
			else
			{
				keepFaces.insert(f);
				allKeepFaces.insert(f);

				// DEBUG:
				//mesh->greenFaces.back().insert(f);
			}

			modifiedFaces.insert(f);
		}

		// Move ownership of cut faces to new vertex
		prev_vd->detachFaces(removeFaces);
		new_vd->attachFaces(removeFaces);
	}

	// recompute normals and reconnect vertex pointers in faces
	mesh->refreshFaces(modifiedFaces);

	mesh->vbo = new VBO(&mesh->vertex, &mesh->vNormal, &mesh->vColor, &mesh->face);
	mesh->setDirtyVBO(true);

	mesh->getUmbrellas();

	return SliceResult(modifiedFaces, cutPointsTable, newPoints, cutEdges, allRemoveFaces, allKeepFaces);
}
