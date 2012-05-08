#pragma once
// Adapted from OpenMesh
#include "GraphicsLibrary/Mesh/QSurfaceMesh.h"

#define get_point(x) (points[x])

/** %Uniform LongestEdgeSubdivision subdivision algorithm
*
* Very simple algorithm splitting all edges which are longer than given via
* set_max_edge_length(). The split is always performed on the longest
* edge in the mesh.
*/
class LongestEdgeSubdivision
{
public:

	typedef double									real_t;
	typedef Surface_mesh							mesh_t;

	typedef std::vector< std::vector<real_t> >      weights_t;
	typedef std::vector<real_t>                     weight_t;

	typedef std::pair<  mesh_t::Edge, real_t >		queueElement;

	class CompareLengthFunction {
	public:
		typedef std::pair< Surface_mesh::Edge, real_t> queueElement;
		bool operator()(const queueElement& t1, const queueElement& t2) // Returns true if t1 is smaller than t2
		{ return (t1.second < t2.second); }
	};

public:

	LongestEdgeSubdivision(double max_edge_len)
	{ set_max_edge_length(max_edge_len); }

public:
	const char *name() const { return "Longest Edge Split"; }

	void set_max_edge_length(double _value) {
		max_edge_length_squared_ = _value * _value;
	}

public:

	bool subdivide( mesh_t& _m, size_t _n , const bool _update_points = true)
	{
		points = _m.vertex_property<Point>("v:point");

		// Sorted queue containing all edges sorted by their decreasing length
		std::priority_queue< queueElement, std::vector< queueElement > , CompareLengthFunction > queue;

		// Build the initial queue
		// First element should be longest edge
		mesh_t::Edge_iterator edgesEnd		= _m.edges_end();
		for (  mesh_t::Edge_iterator eit	= _m.edges_begin(); eit != edgesEnd; ++eit) {
			const  Point to			= points[_m.to_vertex(_m.halfedge(eit,0))];
			const  Point from		= points[_m.from_vertex(_m.halfedge(eit,0))];

			real_t length = (to - from).sqrnorm();

			// Only push the edges that need to be split
			if ( length > max_edge_length_squared_ )
				queue.push( queueElement(eit,length) );
		}

		bool stop = false;
		while ( !stop && ! queue.empty() ) {
			queueElement a = queue.top();
			queue.pop();

			if ( a.second < max_edge_length_squared_ ) {
				stop = true;
				break;
			} else {
				const  Point to   = points[_m.to_vertex(_m.halfedge(a.first,0))];
				const  Point from = points[_m.from_vertex(_m.halfedge(a.first,0))];
				const  Point midpoint = 0.5 * ( to + from );

				const  mesh_t::Vertex newVertex = _m.add_vertex(midpoint);
				_m.split(a.first,newVertex);

				mesh_t::Halfedge_around_vertex_circulator voh_end = _m.halfedges(newVertex);

				for (  mesh_t::Halfedge_around_vertex_circulator voh_it = _m.halfedges(newVertex); voh_it; ) 
				{
					mesh_t::Edge eh		= _m.edge(voh_it);
					const  Point to		= points[_m.to_vertex(voh_it)];
					const  Point from	= points[_m.from_vertex(voh_it)];
					real_t length		= (to - from).sqrnorm();

					// Only push the edges that need to be split
					if ( length > max_edge_length_squared_ )
						queue.push( queueElement(eh,length) );

					++voh_it;
					if(voh_it == voh_end) break;
				}
			}
		}

		return true;
	}

private: // data
	real_t max_edge_length_squared_;

	Surface_mesh::Vertex_property<Point> points;
};
