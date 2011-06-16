#include "Umbrella.h"

Umbrella::Umbrella(VertexDetail * vd)
{
	index = vd->index;
	flag = vd->flag;
	ifaces = vd->ifaces;

	loadHalfEdgeSet();
}

void Umbrella::loadHalfEdgeSet()
{
	halfEdge.clear();

	Face *f, *g;

	// For all incident faces
	for(int i = 0; i < (int)ifaces.size(); i++)
	{
		// Set current 'iface'
		f = ifaces[i];

		// Get two edges connected with us in the current 'iface'
		EdgeSet edges = f->edgesWithVertex(index);

		for(EdgeSet::iterator it = edges.begin(); it != edges.end(); it++)
		{
			int otherVertex = it->neighbor();
			Edge currEdge (index, otherVertex);

			// Second face
			g = NULL;

			// Check if we share an edge by checking all 'ifaces' 
			for(int j = 0; j < (int)ifaces.size(); j++)
			{
				if(i != j && ifaces[j]->hasEdge(currEdge))
				{
					g = ifaces[j];
					break;
				}
			}

			// Insert new half edge into set of half edges
			halfEdge.insert( HalfEdge ( currEdge, f, g ) );
		}
	}

	neighbor.clear();
	neighbor.reserve(halfEdge.size());

	for(HalfEdgeSet::iterator h =  halfEdge.begin(); h != halfEdge.end(); h++)
		neighbor.push_back( (*h).edge.neighbor() );


	// Sort neighbors:
	Vector<int> sorted_neighbor;

	// Faces to keep track of
	Face *start = NULL, *end = NULL, *curr = NULL;

	// Pick a first edge
	HalfEdge he = *halfEdge.begin();
	Edge currEdge(he.edge);

	// First get border check
	PairInt border = borderNeighbours();

	// Having (-1) means no borders
	if(border.first == -1)
	{
		start = he.face(0);
		end = he.face(1);
	
		// If not in order, use the other direction
		if( !he.face(0)->inOrder(he.edge) ){
			start = he.face(1);
			end = he.face(0);
		}
	}
	else
	{
		PairFaces borderFace = borderFaces();

		for(HalfEdgeSet::iterator hi = halfEdge.begin(); hi != halfEdge.end(); hi++)
		{
			// For a border Half Edge
			if(hi->hasFace(NULL))
			{
				if( hi->hasFace(borderFace.first) ) curr = borderFace.first;
				if( hi->hasFace(borderFace.second) ) curr = borderFace.second;

				if(curr && curr->inOrder(hi->edge))
				{
					end = (curr == borderFace.first) ? borderFace.second : borderFace.first;

					// Add border face vertex
					sorted_neighbor.push_back( hi->edge.neighbor() );

					// Move to next face
					PairFaces nei = faceNeighboursPair(curr);
					start = (nei.first != NULL) ? nei.first : nei.second;

					// Set current edge
					currEdge = Edge(index, start->oppositeVertex(Edge(index, sorted_neighbor.back())) );

					sorted_neighbor.push_back(currEdge.neighbor());

					break;
				}
			}
		}
	}

	curr = start;

	while(curr && curr != end)
	{
		sorted_neighbor.push_back( curr->oppositeVertex(currEdge) );

		// Find neighbors of current face
		PairFaces nei = faceNeighboursPair(curr);

		// A good neighbor is one in the correct direction
		Face * nextNeighbour = nei.first;

		if(nei.second && nei.second->hasVertex(sorted_neighbor.back()))
			nextNeighbour = nei.second;
	
		curr = nextNeighbour;
		currEdge = Edge(index, sorted_neighbor.back());
	}

	if(curr)
		sorted_neighbor.push_back( curr->oppositeVertex(currEdge) );

	// replace?
	neighbor = sorted_neighbor;
}

PairInt Umbrella::borderNeighbours()
{
	int result[] = {-1, -1};

	int c = 0;

	for(HalfEdgeSet::iterator h = halfEdge.begin(); h != halfEdge.end(); h++)
	{
		if(h->isBorder())
		{
			result[c++] = h->edge.neighbor();

			if(c > 1)	break;
		}
	}

	if(result[1] == -1) // Non-manifold !!
		result[1] = result[0];

	return PairInt(result[0], result[1]);
}

PairFaces Umbrella::borderFaces()
{
	PairFaces result (NULL, NULL);
	Vector<Face*> ne;

	for(HalfEdgeSet::iterator h = halfEdge.begin(); h != halfEdge.end(); h++)
	{
		if(h->isBorder())
			ne.push_back( h->otherFace(NULL) );
	}

	if(ne.size() > 0) result.first = ne.front();
	if(ne.size() > 1) result.second = ne.back();

	return result;
}

PairFaces Umbrella::sharedFaces(Umbrella * other)
{
	for(HalfEdgeSet::iterator h = halfEdge.begin(); h != halfEdge.end(); h++)
	{
		if(h->edge.neighbor() == other->index)
			return Pair<Face *, Face *>(h->face(0), h->face(1));
	}

	Face * nullFace = 0;

	return PairFaces (nullFace, nullFace);
}

HalfEdge Umbrella::sharedEdge(Umbrella * other)
{
	for(HalfEdgeSet::iterator h = halfEdge.begin(); h != halfEdge.end(); h++)
	{
		if(h->edge.neighbor() == other->index)
			return HalfEdge(*h);
	}

	return HalfEdge();
}

bool Umbrella::hasNeighbour(int index)
{
	for(int i = 0; i < (int)neighbor.size(); i++)
	{
		if(neighbor[i] == index)
			return true;
	}

	return false;
}

bool Umbrella::hasNeighbour(Umbrella * other)
{
	return hasNeighbour(other->index);
}

bool Umbrella::shareTriangle(Umbrella * b, Umbrella * c)
{
	HalfEdge h1 = sharedEdge(b);
	HalfEdge h2 = sharedEdge(c);

	if(h1.isNull() || h2.isNull())
	{
		return false;
	}
	else
	{
		return !(h1.isBorder() && h2.isBorder());
	}
}

int Umbrella::valence()
{
	return neighbor.size();
}

StdSet<int> Umbrella::faceNeighbours(Face * face)
{
	StdSet<int> nighbours;

	for(HalfEdgeSet::iterator h = halfEdge.begin(); h != halfEdge.end(); h++)
	{
		if(h->pFace[0] == face)	nighbours.insert(h->pFace[1]->index);
		if(h->pFace[1] == face)	nighbours.insert(h->pFace[0]->index);
	}

	return nighbours;
}

PairFaces Umbrella::faceNeighboursPair( Face * face )
{
	PairFaces result (NULL, NULL);
	Vector<Face*> ne;

	for(HalfEdgeSet::iterator h = halfEdge.begin(); h != halfEdge.end(); h++)
	{
		if(h->pFace[0] == face)	ne.push_back(h->pFace[1]);
		if(h->pFace[1] == face)	ne.push_back(h->pFace[0]);
	}

	if(ne.size() > 0) result.first = ne.front();
	if(ne.size() > 1) result.second = ne.back();

	return result;
}

Normal Umbrella::normal()
{
	Vec n;
	for(Vector<Face *>::iterator it = ifaces.begin(); it != ifaces.end(); it++)
		n += (*it)->normal();
	n.normalize();
	return n;
}
