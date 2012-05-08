#pragma once
// Adapted from OpenMesh
#include "GraphicsLibrary/Mesh/QSurfaceMesh.h"

class LoopSubdivider
{
public:

	typedef double									real_t;
	typedef Surface_mesh							mesh_t;

	typedef std::pair< real_t, real_t >             weight_t;
	typedef std::vector< std::pair<real_t,real_t> > weights_t;

public:

	LoopSubdivider() : _1over8( 1.0/8.0 ), _3over8( 3.0/8.0 )
	{ init_weights(); }

	~LoopSubdivider() {}

public:

	const char *name() const { return "Uniform Loop"; }

	/// pre-compute weights
	void init_weights(size_t _max_valence=50)
	{
		weights_.resize(_max_valence);
		std::generate(weights_.begin(), weights_.end(), compute_weight());
	}


protected:

	bool prepare( mesh_t& _m )
	{
		points = _m.vertex_property<Point>("v:point");

		vp_pos_ = _m.vertex_property<Point>("v:point2");
		ep_pos_ = _m.edge_property<Point>("e:point");

		return true;
	}

	bool cleanup( mesh_t& _m )
	{
		_m.remove_vertex_property(vp_pos_);
		_m.remove_edge_property(ep_pos_);

		return true;
	}

public:
	bool subdivide( mesh_t& _m, size_t _n, const bool _update_points = true)
	{
		prepare(_m);

		///TODO:Implement fixed positions
		Surface_mesh::Face_iterator   fit, f_end;
		Surface_mesh::Edge_iterator   eit, e_end;
		Surface_mesh::Vertex_iterator vit;

		// Do _n subdivisions
		for (size_t i=0; i < _n; ++i)
		{

			if(_update_points) {
				// compute new positions for old vertices
				for (vit = _m.vertices_begin(); vit != _m.vertices_end(); ++vit) {
					smooth(_m, vit);
				}
			}

			// Compute position for new vertices and store them in the edge property
			for (eit=_m.edges_begin(); eit != _m.edges_end(); ++eit)
				compute_midpoint( _m, eit );

			// Split each edge at midpoint and store precomputed positions (stored in
			// edge property ep_pos_) in the vertex property vp_pos_;

			// Attention! Creating new edges, hence make sure the loop ends correctly.
			e_end = _m.edges_end();
			for (eit=_m.edges_begin(); eit != e_end; ++eit)
				split_edge(_m, eit );


			// Commit changes in topology and reconstitute consistency

			// Attention! Creating new faces, hence make sure the loop ends correctly.
			f_end   = _m.faces_end();
			for (fit = _m.faces_begin(); fit != f_end; ++fit)
				split_face(_m, fit );

			if(_update_points) {
				// Commit changes in geometry
				for ( vit  = _m.vertices_begin(); vit != _m.vertices_end(); ++vit) 
				{
					points[vit] = vp_pos_[vit];
				}
			}
		}

		return cleanup(_m);
	}

private:

	/// Helper functor to compute weights for Loop-subdivision
	/// \internal
	struct compute_weight
	{
		compute_weight() : valence(-1) { }
		weight_t operator() (void)
		{
			//              1
			// alpha(n) = ---- * (40 - ( 3 + 2 cos( 2 Pi / n ) )? )
			//             64

			if (++valence)
			{
				double   inv_v  = 1.0/double(valence);
				double   t      = (3.0 + 2.0 * cos( 2.0 * M_PI * inv_v) );
				double   alpha  = (40.0 - t * t)/64.0;

				return weight_t( 1.0-alpha, inv_v*alpha);
			}
			return weight_t(0.0, 0.0);
		}
		int valence;
	};

private: // topological modifiers

	void split_face(mesh_t& _m, const  Surface_mesh::Face& _fh)
	{
		Surface_mesh::Halfedge
			heh1(_m.halfedge(_fh)),
			heh2(_m.next_halfedge(_m.next_halfedge(heh1))),
			heh3(_m.next_halfedge(_m.next_halfedge(heh2)));

		// Cutting off every corner of the 6_gon
		corner_cutting( _m, heh1 );
		corner_cutting( _m, heh2 );
		corner_cutting( _m, heh3 );
	}


	void corner_cutting(mesh_t& _m, const  Surface_mesh::Halfedge& _he)
	{
		// Define Halfedge Handles
		Surface_mesh::Halfedge
			heh1(_he),
			heh5(heh1),
			heh6(_m.next_halfedge(heh1));

		// Cycle around the polygon to find correct Halfedge
		for (; _m.next_halfedge(_m.next_halfedge(heh5)) != heh1;
			heh5 = _m.next_halfedge(heh5)){}

		Surface_mesh::Vertex
			vh1 = _m.to_vertex(heh1),
			vh2 = _m.to_vertex(heh5);

		Surface_mesh::Halfedge
			heh2(_m.next_halfedge(heh5)),
			heh3(_m.new_edge( vh1, vh2)),
			heh4(_m.opposite_halfedge(heh3));

		/* Intermediate result
		*
		*            *
		*         5 /|\
		*          /_  \
		*    vh2> *     *
		*        /|\3   |\
		*       /_  \|4   \
		*      *----\*----\*
		*          1 ^   6
		*            vh1 (adjust_outgoing halfedge!)
		*/

		// Old and new Face
		Surface_mesh::Face     fh_old(_m.face(heh6));
		Surface_mesh::Face     fh_new(_m.new_face());


		// Re-Set Handles around old Face
		_m.set_next_halfedge(heh4, heh6);
		_m.set_next_halfedge(heh5, heh4);

		_m.set_face(heh4, fh_old);
		_m.set_face(heh5, fh_old);
		_m.set_face(heh6, fh_old);
		_m.set_halfedge(fh_old, heh4);

		// Re-Set Handles around new Face
		_m.set_next_halfedge(heh1, heh3);
		_m.set_next_halfedge(heh3, heh2);

		_m.set_face(heh1, fh_new);
		_m.set_face(heh2, fh_new);
		_m.set_face(heh3, fh_new);

		_m.set_halfedge(fh_new, heh1);
	}


	void split_edge(mesh_t& _m, const  Surface_mesh::Edge& _eh)
	{
		Surface_mesh::Halfedge heh = _m.halfedge(_eh, 0), opp_heh = _m.halfedge(_eh, 1);

		Surface_mesh::Halfedge new_heh, opp_new_heh, t_heh;
		Surface_mesh::Vertex   vh;
		Surface_mesh::Vertex   vh1(_m.to_vertex(heh));
		Point					midP(points[_m.to_vertex(heh)]);
		midP += points[_m.to_vertex(opp_heh)];
		midP *= 0.5;

		// new vertex
		vh = _m.add_vertex(midP);

		// memorize position, will be set later
		vp_pos_[vh] = ep_pos_[_eh];

		// Re-link mesh entities
		if (_m.is_boundary(_eh))
		{
			for (t_heh = heh;
				_m.next_halfedge(t_heh) != opp_heh;
				t_heh = _m.opposite_halfedge(_m.next_halfedge(t_heh))){}
		}
		else
		{
			for (t_heh = _m.next_halfedge(opp_heh);
				_m.next_halfedge(t_heh) != opp_heh;
				t_heh = _m.next_halfedge(t_heh) ){}
		}

		new_heh     = _m.new_edge(vh, vh1);
		opp_new_heh = _m.opposite_halfedge(new_heh);
		_m.set_vertex( heh, vh );

		_m.set_next_halfedge(t_heh, opp_new_heh);
		_m.set_next_halfedge(new_heh, _m.next_halfedge(heh));
		_m.set_next_halfedge(heh, new_heh);
		_m.set_next_halfedge(opp_new_heh, opp_heh);

		if (_m.face(opp_heh).is_valid())
		{
			_m.set_face(opp_new_heh, _m.face(opp_heh));
			_m.set_halfedge(_m.face(opp_new_heh), opp_new_heh);
		}

		_m.set_face( new_heh, _m.face(heh) );
		_m.set_halfedge( vh, new_heh);
		_m.set_halfedge( _m.face(heh), heh );
		_m.set_halfedge( vh1, opp_new_heh );

		// Never forget this, when playing with the topology
		_m.adjust_outgoing_halfedge( vh );
		_m.adjust_outgoing_halfedge( vh1 );
	}

private: // geometry helper

	void compute_midpoint(mesh_t& _m, const  Surface_mesh::Edge& _eh)
	{
		Surface_mesh::Halfedge heh, opp_heh;

		heh      = _m.halfedge( _eh, 0);
		opp_heh  = _m.halfedge( _eh, 1);

		Point pos(points[_m.to_vertex(heh)]);

		pos += points[_m.to_vertex(opp_heh)];

		// boundary edge: just average vertex positions
		if (_m.is_boundary(_eh) )
		{
			pos *= 0.5;
		}
		else // inner edge: add neighboring Vertices to sum
		{
			pos *= real_t(3.0);
			pos += points[_m.to_vertex(_m.next_halfedge(heh))];
			pos += points[_m.to_vertex(_m.next_halfedge(opp_heh))];
			pos *= _1over8;
		}
		ep_pos_[_eh] = pos;
	}

	void smooth(mesh_t& _m, const  Surface_mesh::Vertex& _vh)
	{
		Point pos(0.0,0.0,0.0);

		if (_m.is_boundary(_vh) ) // if boundary: Point 1-6-1
		{
			Surface_mesh::Halfedge heh, prev_heh;
			heh      = _m.halfedge( _vh );

			if ( heh.is_valid() )
			{
				assert( _m.is_boundary( _m.edge( heh ) ) );

				prev_heh = _m.prev_halfedge( heh );

				Surface_mesh::Vertex to_vh = _m.to_vertex( heh ), from_vh = _m.from_vertex( prev_heh );

				// ( v_l + 6 v + v_r ) / 8
				pos  = points[ _vh ];
				pos *= real_t(6.0);
				pos += points[ to_vh ];
				pos += points[ from_vh ];
				pos *= _1over8;

			}
			else
				return;
		}
		else // inner vertex: (1-a) * p + a/n * Sum q, q in one-ring of p
		{
			size_t valence(0);

			Surface_mesh::Vertex_around_vertex_circulator vit, vend;
			vit = vend = _m.vertices(_vh);
			do{ 
				++valence;
				pos += points[vit];
				++vit;
			} while(vit != vend);

			pos *= weights_[valence].second; // alpha(n)/n * Sum q, q in one-ring of p
			pos += weights_[valence].first * points[_vh]; // + (1-a)*p
		}

		vp_pos_[_vh] = pos;
	}

private: // data

	Surface_mesh::Vertex_property<Point> vp_pos_;
	Surface_mesh::Edge_property<Point> ep_pos_;

	Surface_mesh::Vertex_property<Point> points;

	weights_t     weights_;

	const real_t _1over8;
	const real_t _3over8;
};
