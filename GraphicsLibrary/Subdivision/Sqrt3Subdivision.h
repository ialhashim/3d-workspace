#pragma once
// Adapted from OpenMesh
#include "GraphicsLibrary/Mesh/QSurfaceMesh.h"

#if 0
#define get_point(x) (points[x])
#define set_point(x) points[x]=


class Sqrt3T
{
public:
	typedef double									real_t;
	typedef QSurfaceMesh							mesh_t;

	typedef std::pair< real_t, real_t >             weight_t;
	typedef std::vector< std::pair<real_t,real_t> > weights_t;

public:

	Sqrt3T(mesh_t &_m): _1over3( 1.0/3.0 ), _1over27( 1.0/27.0 )
	{ init_weights(); }

public:

	const char *name() const { return "Uniform Sqrt3"; }

	/// Pre-compute weights
	void init_weights(size_t _max_valence=50)
	{
		weights_.resize(_max_valence);
		std::generate(weights_.begin(), weights_.end(), compute_weight());
	}

protected:

	bool prepare( Surface_mesh& _m )
	{
		points = _m.vertex_property<Point>("v:point");

		vp_pos_ = _m.vertex_property<Point>("v:point2");
		ep_nv_ = _m.edge_property< std::pair< Surface_mesh::Vertex,	 Surface_mesh::Vertex> >("e:nv");
		mp_gen_ = 0;

		return true;
	}
	
	bool cleanup( Surface_mesh& _m )
	{
		_m.remove_vertex_property(vp_pos_);
		_m.remove_edge_property( ep_nv_ );
		return true;
	}

public:
	bool subdivide( Surface_mesh& _m, size_t _n , const bool _update_points = true)
	{
		prepare(_m);

		///TODO:Implement fixed positions
		Surface_mesh::Vertex_iterator			vit;
		Surface_mesh::VertexVertex_iterator		vvit;
		Surface_mesh::Edge_iterator				eit;
		Surface_mesh::Face_iterator				fit;
		Surface_mesh::FaceVertex_iterator		fvit;
		Surface_mesh::Vertex					vh;
		Surface_mesh::Halfedge					heh;
		Surface_mesh::Point						pos(0,0,0), zero(0,0,0);
		size_t									&gen = _m.property( mp_gen_ );

		for (size_t l=0; l<_n; ++l)
		{
			// tag existing edges
			for (eit=_m.edges_begin(); eit != _m.edges_end();++eit)
			{
				_m.status( eit ).set_tagged( true );
				if ( (gen%2) && _m.is_boundary(eit) )
					compute_new_boundary_points( _m, eit ); // *) creates new vertices
			}

			// do relaxation of old vertices, but store new pos in property vp_pos_

			for (vit=_m.vertices_begin(); vit!=_m.vertices_end(); ++vit)
			{
				if ( _m.is_boundary(vit) )
				{
					if ( gen%2 )
					{
						heh  = _m.halfedge(vit);
						if (heh.is_valid()) // skip isolated newly inserted vertices *)
						{
							Surface_mesh::Halfedge prev_heh = _m.prev_halfedge(heh);

							assert( _m.is_boundary(heh     ) );
							assert( _m.is_boundary(prev_heh) );

							pos  = get_point(_m.to_vertex(heh));
							pos += get_point(_m.from_vertex(prev_heh));
							pos *= real_t(4.0);

							pos += real_t(19.0) * get_point( vit );
							pos *= _1over27;

							_m.property( vp_pos_, vit ) = pos;
						}
					}
					else
						_m.property( vp_pos_, vit ) = get_point( vit );
				}
				else
				{
					size_t valence=0;

					pos = zero;
					for ( vvit = _m.vv_iter(vit); vvit; ++vvit)
					{
						pos += get_point( vvit );
						++valence;
					}
					pos *= weights_[ valence ].second;
					pos += weights_[ valence ].first * get_point(vit);
					_m.property( vp_pos_, vit ) =  pos;
				}
			}   

			// insert new vertices, but store pos in vp_pos_
			Surface_mesh::Face_iterator fend = _m.faces_end();
			for (fit = _m.faces_begin();fit != fend; ++fit)
			{
				if ( (gen%2) && _m.is_boundary(fit))
				{
					boundary_split( _m, fit );
				}
				else
				{
					fvit = _m.fv_iter( fit );        
					pos  = get_point(  fvit);
					pos += get_point(++fvit);
					pos += get_point(++fvit);
					pos *= _1over3;
					vh   = _m.add_vertex( zero );
					_m.property( vp_pos_, vh ) = pos;
					_m.split( fit, vh );
				}
			}

			// commit new positions (now iterating over all vertices)
			for (vit=_m.vertices_begin();vit != _m.vertices_end(); ++vit)
				set_point(vit, _m.property( vp_pos_, vit ) );

			// flip old edges
			for (eit=_m.edges_begin(); eit != _m.edges_end(); ++eit)
				if ( _m.status( eit ).tagged() && !_m.is_boundary( eit ) )
					_m.flip(eit);

			// Now we have an consistent mesh!
			//ASSERT_CONSISTENCY( Surface_mesh, _m );

			// increase generation by one
			++gen;
		}

		return cleanup(_m);
	}

private:

	/// Helper functor to compute weights for sqrt(3)-subdivision
	/// \internal
	struct compute_weight 
	{
		compute_weight() : valence(-1) { }    
		weight_t operator() (void) 
		{ 
			if (++valence)
			{
				real_t alpha = (4.0-2.0*cos(2.0*M_PI / (double)valence))/9.0;
				return weight_t( real_t(1)-alpha, alpha/real_t(valence) );
			}
			return weight_t(0.0, 0.0);
		}    
		int valence;
	};

private:

	// Pre-compute location of new boundary points for odd generations
	// and store them in the edge property ep_nv_;
	void compute_new_boundary_points( Surface_mesh& _m, const  Surface_mesh::EdgeHandle& _eh)
	{
		assert( _m.is_boundary(_eh) );

		Surface_mesh::Halfedge heh;
		Surface_mesh::Vertex   vh1, vh2, vh3, vh4, vhl, vhr;
		Point zero(0,0,0), P1, P2, P3, P4;

		/*
		//       *---------*---------*
		//      / \       / \       / \
		//     /   \     /   \     /   \
		//    /     \   /     \   /     \
		//   /       \ /       \ /       \
		//  *---------*--#---#--*---------*
		//                
		//  ^         ^  ^   ^  ^         ^
		//  P1        P2 pl  pr P3        P4
		*/
		// get halfedge pointing from P3 to P2 (outer boundary halfedge)

		heh = _m.halfedge(_eh, _m.is_boundary(_m.halfedge(_eh,1)));

		assert( _m.is_boundary( _m.next_halfedge( heh ) ) );
		assert( _m.is_boundary( _m.prev_halfedge( heh ) ) );

		vh1 = _m.to_vertex( _m.next_halfedge( heh ) );
		vh2 = _m.to_vertex( heh );
		vh3 = _m.from_vertex( heh );
		vh4 = _m.from_vertex( _m.prev_halfedge( heh ));

		P1  = get_point(vh1);
		P2  = get_point(vh2);
		P3  = get_point(vh3);
		P4  = get_point(vh4);

		vhl = _m.add_vertex(zero);
		vhr = _m.add_vertex(zero);

		_m.property(vp_pos_, vhl ) = (P1 + real_t(16.0f) * P2 + real_t(10.0f) * P3) * _1over27;
		_m.property(vp_pos_, vhr ) = ( real_t(10.0f) * P2 + real_t(16.0f) * P3 + P4) * _1over27;
		_m.property(ep_nv_, _eh).first  = vhl;
		_m.property(ep_nv_, _eh).second = vhr; 
	}


	void boundary_split( Surface_mesh& _m, const  Surface_mesh::FaceHandle& _fh )
	{
		assert( _m.is_boundary(_fh) );

		Surface_mesh::Vertex     vhl, vhr;
		Surface_mesh::FaceEdge_iterator     fe_it;
		Surface_mesh::Halfedge   heh;

		// find boundary edge
		for( fe_it=_m.fe_iter( _fh ); fe_it && !_m.is_boundary( fe_it ); ++fe_it ) {};

		// use precomputed, already inserted but not linked vertices
		vhl = _m.property(ep_nv_, fe_it).first;
		vhr = _m.property(ep_nv_, fe_it).second;

		/*
		//       *---------*---------*
		//      / \       / \       / \
		//     /   \     /   \     /   \
		//    /     \   /     \   /     \
		//   /       \ /       \ /       \
		//  *---------*--#---#--*---------*
		//                
		//  ^         ^  ^   ^  ^         ^
		//  P1        P2 pl  pr P3        P4
		*/
		// get halfedge pointing from P2 to P3 (inner boundary halfedge)

		heh = _m.halfedge(fe_it, _m.is_boundary(_m.halfedge(fe_it,0)));

		Surface_mesh::Halfedge pl_P3;

		// split P2->P3 (heh) in P2->pl (heh) and pl->P3
		boundary_split( _m, heh, vhl );         // split edge
		pl_P3 = _m.next_halfedge( heh ); // store next halfedge handle
		boundary_split( _m, heh );              // split face

		// split pl->P3 in pl->pr and pr->P3
		boundary_split( _m, pl_P3, vhr );
		boundary_split( _m, pl_P3 );

		assert( _m.is_boundary( vhl ) && _m.halfedge(vhl).is_valid() );
		assert( _m.is_boundary( vhr ) && _m.halfedge(vhr).is_valid() );
	}

	void boundary_split(Surface_mesh& _m, const  Surface_mesh::Halfedge& _heh, const  Surface_mesh::Vertex& _vh)
	{
		assert( _m.is_boundary( _m.edge(_heh) ) );

		Surface_mesh::Halfedge heh(_heh), opp_heh( _m.opposite_halfedge(_heh) ), new_heh, opp_new_heh;
		Surface_mesh::Vertex   to_vh(_m.to_vertex(heh));
		Surface_mesh::Halfedge t_heh;

		/*
		*            P5
		*             *
		*            /|\
		*           /   \
		*          /     \
		*         /       \
		*        /         \
		*       /_ heh  new \
		*      *-----\*-----\*\-----*
		*             ^      ^ t_heh
		*            _vh     to_vh
		*
		*     P1     P2     P3     P4
		*/
		// Re-Setting Handles

		// find halfedge point from P4 to P3
		for(t_heh = heh; _m.next_halfedge(t_heh) != opp_heh; t_heh = _m.opposite_halfedge(_m.next_halfedge(t_heh))){}

		assert( _m.is_boundary( t_heh ) );

		new_heh     = _m.new_edge( _vh, to_vh );
		opp_new_heh = _m.opposite_halfedge(new_heh);

		// update halfedge connectivity

		_m.set_next_halfedge(t_heh,   opp_new_heh); // P4-P3 -> P3-P2
		// P2-P3 -> P3-P5
		_m.set_next_halfedge(new_heh, _m.next_halfedge(heh));    
		_m.set_next_halfedge(heh,         new_heh); // P1-P2 -> P2-P3
		_m.set_next_halfedge(opp_new_heh, opp_heh); // P3-P2 -> P2-P1

		// both opposite halfedges point to same face
		_m.set_face(opp_new_heh, _m.face(opp_heh));

		// let heh finally point to new inserted vertex
		_m.set_vertex(heh, _vh); 

		// let heh and new_heh point to same face
		_m.set_face(new_heh, _m.face(heh));

		// let opp_new_heh be the new outgoing halfedge for to_vh 
		// (replaces for opp_heh)
		_m.set_halfedge( to_vh, opp_new_heh );

		// let opp_heh be the outgoing halfedge for _vh
		_m.set_halfedge( _vh, opp_heh );
	}

	void boundary_split( Surface_mesh& _m, const  Surface_mesh::Halfedge& _heh)
	{
		assert( _m.is_boundary( _m.opposite_halfedge( _heh ) ) );

		 Surface_mesh::Halfedge 
			heh(_heh),
			n_heh(_m.next_halfedge(heh));

		 Surface_mesh::Vertex   
			to_vh(_m.to_vertex(heh));

		 Surface_mesh::Halfedge 
			heh2(_m.new_edge(to_vh,	_m.to_vertex(_m.next_halfedge(n_heh)))),
			heh3(_m.opposite_halfedge(heh2));

		 Surface_mesh::Face
			new_fh(_m.new_face()),
			fh(_m.face(heh));

		// Relink (half)edges    
#define set_next_heh set_next_halfedge
#define next_heh next_halfedge

		_m.set_face(heh,  new_fh);
		_m.set_face(heh2, new_fh);
		_m.set_next_heh(heh2, _m.next_heh(_m.next_heh(n_heh)));
		_m.set_next_heh(heh,  heh2);
		_m.set_face( _m.next_heh(heh2), new_fh);

		// _m.set_face( _m.next_heh(_m.next_heh(heh2)), new_fh);

		_m.set_next_heh(heh3,                           n_heh);
		_m.set_next_heh(_m.next_halfedge(n_heh), heh3);
		_m.set_face(heh3,  fh);
		// _m.set_face(n_heh, fh);

		_m.set_halfedge(    fh, n_heh);
		_m.set_halfedge(new_fh, heh);

#undef set_next_halfedge
#undef next_halfedge

	}

private:

	weights_t     weights_;

	Surface_mesh::Vertex_property<Point> vp_pos_;
	Surface_mesh::Edge_property<std::pair< Surface_mesh::Vertex,	 Surface_mesh::Vertex> > ep_nv_;

	Surface_mesh::Vertex_property<Point> points;

	size_t mp_gen_;

	const real_t _1over3;
	const real_t _1over27;
};
#endif
