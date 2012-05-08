#pragma once 
// Adapted from http://code.google.com/p/starlab/

#include "GraphicsLibrary/Mesh/QSurfaceMesh.h"
typedef double Scalar;
typedef int Counter;

class LaplacianRemesher{
private:
    Scalar target_edge_length;
	Surface_mesh * mesh;
	Surface_mesh::Vertex_property<Point> points;
	Surface_mesh::Vertex_property<Normal> vnormal;

public:
    LaplacianRemesher(Surface_mesh* mesh, Scalar target_edge_length)
	{
		this->mesh = mesh;
        this->target_edge_length = target_edge_length;

		points = mesh->vertex_property<Point>("v:point");
		vnormal = mesh->vertex_property<Point>("v:normal");
    }
  
    void tangential_smoothing(int num_iter)
	{
        Surface_mesh::Vertex_iterator v_it, v_end(mesh->vertices_end());
        Surface_mesh::Vertex_around_vertex_circulator  vv_it, vv_end;
        Scalar  valence;
        Point   u;
    
        for (Counter iters=0; iters<num_iter; ++iters){
            // normals have to be ok at the beginning of each iteration!
            
            // compute tangential updates
            for (v_it=mesh->vertices_begin(); v_it!=v_end; ++v_it)
            {
                if (!mesh->is_boundary(v_it))
                {
                    // compute barycenter of neighbors
                    u = 0.0;
                    valence = 0;
                    vv_it = vv_end = mesh->vertices(v_it);
                    do
                    {
                        u += points[vv_it];
                        ++valence;
                    }
                    while (++vv_it != vv_end);
                    u *= (1.0/valence);
                    
                    // project vector from point to barycenter to tangent plane
                    u -= points[v_it];                
                    u -= vnormal[v_it] * dot(u, vnormal[v_it]);
                    
                    // move point
                    points[v_it] += u;
                }
            }
        
            // update normals
            mesh->update_vertex_normals();
        }
    }
    
    void split_long_edges(Counter numiters){
        /** split edges that are too long
         \li loop several times over all edges (note that mesh->edges_end() might be different in each iteration!).
         \li determine whether the edge is too long using is_too_long().
         \li if an edge is too long, split it (see Surface_mesh::split(Edge, Vertex)), 
         set the position of the new vertex to the edge midpoint, and set its 
         normal to the average of the endpoint normals.
         */
        Surface_mesh::Edge_iterator e_it, e_end;
        Surface_mesh::Vertex v, v0, v1;
        bool done;
        Counter i;
        
        // do at most 10 loops over all edges
        for (done=false, i=0; !done && i<numiters; ++i)
        {
            done = true;
            
            for (e_it=mesh->edges_begin(), e_end=mesh->edges_end(); e_it!=e_end; ++e_it)
            {
                v0 = mesh->vertex(e_it, 0);
                v1 = mesh->vertex(e_it, 1);
                
                if(is_too_long(v0, v1))
                {                
                    v = mesh->add_vertex((points[v0] + points[v1]) * 0.5);
                    vnormal[v] = (vnormal[v0] + vnormal[v1]).normalize();
                    mesh->split(e_it, v);
                    done = false;
                }
            }
        }
    }
    
    void collapse_short_edges(Counter numiters){
        Surface_mesh::Edge_iterator     eit, eend;
        Surface_mesh::Vertex_around_vertex_circulator vv_it, vv_end;
        Surface_mesh::Vertex    v0, v1;
        Surface_mesh::Halfedge  h;
        bool  done;
        Counter   i;
        
        for (done=false, i=0; !done && i<numiters; ++i){
            done = true;
            
            for (eit=mesh->edges_begin(), eend=mesh->edges_end(); eit!=eend; ++eit)
            {
                h   = mesh->halfedge(eit, 0);
                v0  = mesh->from_vertex(h);
                v1  = mesh->to_vertex(h);
                    
    
                // edge too short -> we want to collapse
                if (is_too_short(v0, v1))
                {
                    // is collapse allowed?
                    if (mesh->is_collapse_ok(h))
                    {                    
                        // check that no too long edge will be created
                        bool ok=true;
                        vv_it = vv_end = mesh->vertices(v0);
                        do 
                        {
                            if (is_too_long(vv_it, v1))
                            {
                                ok = false;
                                break;
                            }
                        } 
                        while (++vv_it != vv_end);
                        
                        if (ok)
                        {
                            mesh->collapse(h);                    
                            done = false;
                        }
                    }
                }
            }
        }
        
        mesh->garbage_collection();
    }

    /// @brief flip edges to improve vertex valences
    void flip_edges(Counter num_edge_flips){
        /** 
         \li loop several times over all edges (only consider non-boundary edges)
         \li collect the four vertices of the two incident triangles
         \li determine the optimal valences for these vertices (boundary or non-boundary vertex?)
         \li if an edge flip improves the valences, then flip it
         \li compare the sum of squared deviations from the optimal values 
         */
        
        Surface_mesh::Edge_iterator     eit, eend;
        Surface_mesh::Vertex            v0, v1, v2, v3;
        Surface_mesh::Halfedge          h;
        int                             val0, val1, val2, val3;
        int                             val_opt0, val_opt1, val_opt2, val_opt3;
        int                             ve_before, ve_after;
        bool                            done;
        Counter                             i;
        
        for (done=false, i=0; !done && i<num_edge_flips; ++i){
            done = true;
            
            for (eit=mesh->edges_begin(), eend=mesh->edges_end(); eit!=eend; ++eit){
                if (!mesh->is_boundary(eit)){
                    h  = mesh->halfedge(eit, 0);
                    v0 = mesh->to_vertex(h);
                    v2 = mesh->to_vertex(mesh->next_halfedge(h));
                    h  = mesh->halfedge(eit, 1);
                    v1 = mesh->to_vertex(h);
                    v3 = mesh->to_vertex(mesh->next_halfedge(h));
                    
                    val0 = mesh->valence(v0);
                    val1 = mesh->valence(v1);
                    val2 = mesh->valence(v2);
                    val3 = mesh->valence(v3);
                    
                    val_opt0 = (mesh->is_boundary(v0) ? 4 : 6);
                    val_opt1 = (mesh->is_boundary(v1) ? 4 : 6);
                    val_opt2 = (mesh->is_boundary(v2) ? 4 : 6);
                    val_opt3 = (mesh->is_boundary(v3) ? 4 : 6);
                    
                    ve_before = (pow(val0 - val_opt0, 2.0f) +
                                 pow(val1 - val_opt1, 2.0f) +
                                 pow(val2 - val_opt2, 2.0f) +
                                 pow(val3 - val_opt3, 2.0f));
                    
                    --val0;  --val1;
                    ++val2;  ++val3;
                    
                    ve_after = (pow(val0 - val_opt0, 2.0f) +
                                pow(val1 - val_opt1, 2.0f) +
                                pow(val2 - val_opt2, 2.0f) +
                                pow(val3 - val_opt3, 2.0f));
                    
                    if (ve_before > ve_after && mesh->is_flip_ok(eit)){
                        mesh->flip(eit);
                        done = false;
                    }
                }
            }
        }
    }
    
private:
    /// returns whether the edge (v0,v1) is too long
    bool is_too_long(Surface_mesh::Vertex v0, Surface_mesh::Vertex v1) const
    { 
        return (points[v0]-points[v1]).norm() > 4.0/3.0*target_edge_length; 
    }

    /// returns whether the edge (v0,v1) is too short
    bool is_too_short(Surface_mesh::Vertex v0, Surface_mesh::Vertex v1) const
    { 
        return (points[v0]-points[v1]).norm() < 4.0/5.0*target_edge_length; 
    }

public:
	static void remesh(QSurfaceMesh *model, double longest_edge_length, int iter = 2, int smoothIter = 3)
	{
		remesh(model, longest_edge_length, iter, iter, iter, smoothIter);
	}

	static void remesh(Surface_mesh *mesh, Scalar longest_edge_length, Counter num_split_iters, 
		Counter num_collapse_iters, Counter num_edge_flips, Counter num_smooth_iters)
	{
		LaplacianRemesher h(mesh, longest_edge_length);

		/// Perform refinement
		h.split_long_edges(num_split_iters);
		h.collapse_short_edges(num_collapse_iters);
		h.flip_edges(num_edge_flips);
		h.tangential_smoothing(num_smooth_iters);
	}
};
