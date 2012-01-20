#include "Curvature.h"

// Rotate a coordinate system to be perpendicular to the given normal
void Curvature::rot_coord_sys( const Point &old_u, const Point &old_v, const Point &new_norm, Point &new_u, Point &new_v )
{
	new_u = old_u;
	new_v = old_v;
	Point old_norm = cross(old_u, old_v);
	double ndot = dot(old_norm, new_norm);
	if (unlikely(ndot <= -1.0f)) {
		new_u = -new_u;
		new_v = -new_v;
		return;
	}
	Point perp_old = new_norm - ndot * old_norm;
	Point dperp = 1.0f / (1 + ndot) * (old_norm + new_norm);
	new_u -= dperp * dot(new_u , perp_old);
	new_v -= dperp * dot(new_v , perp_old);
}

// Re-project a curvature tensor from the basis spanned by old_u and old_v
// (which are assumed to be unit-length and perpendicular) to the
// new_u, new_v basis.
void Curvature::proj_curv( const Point &old_u, const Point &old_v, double old_ku, double old_kuv, double old_kv, 
	const Point &new_u, const Point &new_v, double &new_ku, double &new_kuv, double &new_kv )
{
	Point r_new_u, r_new_v;
	rot_coord_sys(new_u, new_v, cross(old_u, old_v), r_new_u, r_new_v);

	double u1 = dot(r_new_u, old_u);
	double v1 = dot(r_new_u, old_v);
	double u2 = dot(r_new_v, old_u);
	double v2 = dot(r_new_v, old_v);

	new_ku = old_ku * u1*u1 + old_kuv * (2.0f * u1*v1) + old_kv * v1*v1;
	new_kuv = old_ku * u1*u2 + old_kuv * (u1*v2 + u2*v1) + old_kv * v1*v2;
	new_kv = old_ku * u2*u2 + old_kuv * (2.0f * u2*v2) + old_kv * v2*v2;
}

// And for derivatives of curvature
void Curvature::proj_dcurv(const Point &old_u, const Point &old_v, const Vec4d old_dcurv, 
	const Point &new_u, const Point &new_v, Vec4d & new_dcurv)
{
	Point r_new_u, r_new_v;
	rot_coord_sys(new_u, new_v, cross(old_u, old_v), r_new_u, r_new_v);

	double u1 = dot(r_new_u, old_u);
	double v1 = dot(r_new_u, old_v);
	double u2 = dot(r_new_v, old_u);
	double v2 = dot(r_new_v, old_v);

	new_dcurv[0] = old_dcurv[0]*u1*u1*u1 + old_dcurv[1]*3.0f*u1*u1*v1 + old_dcurv[2]*3.0f*u1*v1*v1 + old_dcurv[3]*v1*v1*v1;
	new_dcurv[1] = old_dcurv[0]*u1*u1*u2 + old_dcurv[1]*(u1*u1*v2 + 2.0f*u2*u1*v1) + old_dcurv[2]*(u2*v1*v1 + 2.0f*u1*v1*v2) + old_dcurv[3]*v1*v1*v2;
	new_dcurv[2] = old_dcurv[0]*u1*u2*u2 + old_dcurv[1]*(u2*u2*v1 + 2.0f*u1*u2*v2) + old_dcurv[2]*(u1*v2*v2 + 2.0f*u2*v2*v1) + old_dcurv[3]*v1*v2*v2;
	new_dcurv[3] = old_dcurv[0]*u2*u2*u2 + old_dcurv[1]*3.0f*u2*u2*v2 + old_dcurv[2]*3.0f*u2*v2*v2 + old_dcurv[3]*v2*v2*v2;
}

// Given a curvature tensor, find principal directions and curvatures
// Makes sure that pdir1 and pdir2 are perpendicular to normal
void Curvature::diagonalize_curv( const Point &old_u, const Point &old_v, double ku, double kuv, double kv, 
	const Point &new_norm, Point &pdir1, Point &pdir2, double &k1, double &k2 )
{
	Point r_old_u, r_old_v;
	rot_coord_sys(old_u, old_v, new_norm, r_old_u, r_old_v);

	double c = 1, s = 0, tt = 0;

	if (likely(kuv != 0.0f)) 
	{
		// Jacobi rotation to diagonalize
		double h = 0.5f * (kv - ku) / kuv;
		tt = (h < 0.0f) ? 1.0f / (h - sqrt(1.0f + h*h)) : 1.0f / (h + sqrt(1.0f + h*h));
		c = 1.0f / sqrt(1.0f + tt*tt);
		s = tt * c;
	}

	k1 = ku - tt * kuv;
	k2 = kv + tt * kuv;

	if (fabs(k1) >= fabs(k2)) {
		pdir1 = c*r_old_u - s*r_old_v;
	} else {
		SWAP(k1, k2, double);
		pdir1 = s*r_old_u + c*r_old_v;
	}

	pdir2 = cross(new_norm, pdir1);
}

void Curvature::computePrincipalCurvatures( Surface_mesh * src_mesh )
{
	computePointAreas(src_mesh);

	uint nv = src_mesh->n_vertices(), nf = src_mesh->n_faces();

	curv1.resize(nv);
	curv2.resize(nv);

	pdir1.resize(nv);
	pdir2.resize(nv);

	std::vector<double> curv12(nv);
	
	Surface_mesh::Vertex_property<Point> points = src_mesh->vertex_property<Point>("v:point");
	Surface_mesh::Vertex_property<Normal> normals = src_mesh->vertex_property<Normal>("v:normal");
	Surface_mesh::Face_iterator fit, fend = src_mesh->faces_end();
	Surface_mesh::Vertex_around_face_circulator fvit, fvend;
	Surface_mesh::Vertex_iterator vit, vend = src_mesh->vertices_end();

	//#pragma omp parallel for
	for (fit = src_mesh->faces_begin(); fit != fend; ++fit) 
	{
		fvit = src_mesh->vertices(fit);

		uint vi0 = ((Vertex)fvit).idx();  Point v0 = points[fvit]; ++fvit;
		uint vi1 = ((Vertex)fvit).idx();  Point v1 = points[fvit]; ++fvit;
		uint vi2 = ((Vertex)fvit).idx();  Point v2 = points[fvit];

		pdir1[vi0] = v1 - v0;
		pdir1[vi1] = v2 - v1;
		pdir1[vi2] = v0 - v2;
	}

	//#pragma omp parallel for
	for (vit = src_mesh->vertices_begin(); vit != vend; ++vit) 
	{
		uint vi = ((Vertex)vit).idx();

		pdir1[vi] = cross(pdir1[vi], normals[vit]);
		pdir1[vi].normalize();
		pdir2[vi] = cross(normals[vit], pdir1[vi]);
	}

	// Compute curvature per-face
	//#pragma omp parallel for
	for (fit = src_mesh->faces_begin(); fit != fend; ++fit) 
	{
		uint i = ((Surface_mesh::Face)fit).idx(); // face index

		uint vi[3];	Normal vn[3];
		fvit = src_mesh->vertices(fit);
		
		vi[0] = ((Vertex)fvit).idx(); Point v0 = points[fvit]; vn[0] = normals[fvit]; ++fvit;
		vi[1]  = ((Vertex)fvit).idx(); Point v1 = points[fvit]; vn[1] = normals[fvit]; ++fvit;
		vi[2]  = ((Vertex)fvit).idx(); Point v2 = points[fvit]; vn[2] = normals[fvit]; ++fvit;

		// Edges
		Point e[3] = {v2 - v1,  v0 - v2,  v1 - v0};

		// N-T-B coordinate system per face
		Point t = e[0];
		t.normalize();
		Point n = cross(e[0], e[1]);
		Point b = cross(n, t);
		b.normalize();

		// Estimate curvature based on variation of normals
		// along edges
		double m[3] = { 0, 0, 0 };
		double w[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };
		for (int j = 0; j < 3; j++) {
			double u = dot(e[j], t);
			double v = dot(e[j], b);
			w[0][0] += u*u;
			w[0][1] += u*v;
			//w[1][1] += v*v + u*u; 
			//w[1][2] += u*v; 
			w[2][2] += v*v;
			Point dn = vn[PREV_Index(j)] - vn[NEXT_Index(j)];
			double dnu = dot(dn, t);
			double dnv = dot(dn, b);
			m[0] += dnu*u;
			m[1] += dnu*v + dnv*u;
			m[2] += dnv*v;
		}
		w[1][1] = w[0][0] + w[2][2];
		w[1][2] = w[0][1];

		// Least squares solution
		double diag[3];
		if (!ldltdc<double,3>(w, diag)) {
			//printf("ldltdc failed!\n");
			continue;
		}
		ldltsl<double,3>(w, diag, m, m);

		// Push it back out to the vertices
		for (uint j = 0; j < 3; j++) 
		{
			int vj = vi[j];
			double c1, c12, c2;
			proj_curv(t, b, m[0], m[1], m[2], pdir1[vj], pdir2[vj], c1, c12, c2);
			double wt = cornerareas[i][j] / pointareas[vj];

			//#pragma omp atomic
			curv1[vj] += wt * c1;

			//#pragma omp atomic
			curv12[vj] += wt * c12;

			//#pragma omp atomic
			curv2[vj] += wt * c2;
		}
	}

	// Store results into Surface_mesh object
	Surface_mesh::Vertex_property<double> my_curv1 = src_mesh->vertex_property<double>("v:curv1");
	Surface_mesh::Vertex_property<double> my_curv2 = src_mesh->vertex_property<double>("v:vurv2");
	Surface_mesh::Vertex_property<Vec3d> my_pdir1 = src_mesh->vertex_property<Vec3d>("v:pdir1");
	Surface_mesh::Vertex_property<Vec3d> my_pdir2 = src_mesh->vertex_property<Vec3d>("v:pdir2");

	// Compute principal directions and curvatures at each vertex
	//#pragma omp parallel for
	for (vit = src_mesh->vertices_begin(); vit != vend; ++vit) 
	{
		uint i = ((Surface_mesh::Vertex)vit).idx();
		diagonalize_curv(pdir1[i], pdir2[i],curv1[i], curv12[i], curv2[i], normals[vit], pdir1[i], pdir2[i],curv1[i], curv2[i]);

		my_curv1[vit] = curv1[i];
		my_curv2[vit] = curv2[i];

		my_pdir1[vit] = pdir1[i];
		my_pdir2[vit] = pdir2[i];
	}
}

// Compute derivatives of curvature
void Curvature::computeDerivativesCurvatures(Surface_mesh * src_mesh)
{
	// Compute principal curvatures and directions
	computePrincipalCurvatures(src_mesh);

	Surface_mesh::Vertex_property<Point> points = src_mesh->vertex_property<Point>("v:point");
	Surface_mesh::Vertex_property<Normal> normals = src_mesh->vertex_property<Normal>("v:normal");
	Surface_mesh::Face_iterator fit, fend = src_mesh->faces_end();
	Surface_mesh::Vertex_around_face_circulator fvit, fvend;

	// Resize the arrays we'll be using
	uint nv = src_mesh->n_vertices(), nf = src_mesh->n_faces();

	dcurv.resize(nv);

	// Compute derivatives of curvature per-face
	//#pragma omp parallel for
	for (fit = src_mesh->faces_begin(); fit != fend; ++fit) 
	{
		uint i = ((Surface_mesh::Face)fit).idx(); // face index

		uint vi[3];	Normal vn[3];
		fvit = src_mesh->vertices(fit);

		vi[0] = ((Vertex)fvit).idx(); Point v0 = points[fvit]; vn[0] = normals[fvit]; ++fvit;
		vi[1]  = ((Vertex)fvit).idx(); Point v1 = points[fvit]; vn[1] = normals[fvit]; ++fvit;
		vi[2]  = ((Vertex)fvit).idx(); Point v2 = points[fvit]; vn[2] = normals[fvit]; ++fvit;

		// Edges
		Point e[3] = {v2 - v1,  v0 - v2,  v1 - v0};

		// N-T-B coordinate system per face
		Point t = e[0];
		t.normalize();
		Point n = cross(e[0], e[1]);
		Point b = cross(n, t);
		b.normalize();

		// Project curvature tensor from each vertex into this
		// face's coordinate system
		Point fcurv[3];
		for (int j = 0; j < 3; j++) 
		{
			int vj = vi[j];
			proj_curv(pdir1[vj], pdir2[vj], curv1[vj], 0, curv2[vj],t, b, fcurv[j][0], fcurv[j][1], fcurv[j][2]);
		}

		// Estimate derivatives of curvature based on variation of curvature along edges
		double m[4] = {0,0,0,0};
		double w[4][4] = { {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0} };

		for (int j = 0; j < 3; j++) 
		{
			// Variation of curvature along each edge
			Point dfcurv = fcurv[PREV_Index(j)] - fcurv[NEXT_Index(j)];
			double u = dot(e[j], t);
			double v = dot(e[j], b);
			double u2 = u*u, v2 = v*v, uv = u*v;
			w[0][0] += u2;
			w[0][1] += uv;
			//w[1][1] += 2.0f*u2 + v2;
			//w[1][2] += 2.0f*uv;
			//w[2][2] += u2 + 2.0f*v2;
			//w[2][3] += uv;
			w[3][3] += v2;
			m[0] += u*dfcurv[0];
			m[1] += v*dfcurv[0] + 2.0f*u*dfcurv[1];
			m[2] += 2.0f*v*dfcurv[1] + u*dfcurv[2];
			m[3] += v*dfcurv[2];
		}

		w[1][1] = 2.0f * w[0][0] + w[3][3];
		w[1][2] = 2.0f * w[0][1];
		w[2][2] = w[0][0] + 2.0f * w[3][3];
		w[2][3] = w[0][1];

		// Least squares solution
		double d[4];

		if (!ldltdc<double,4>(w, d)) 
		{
			//printf("ldltdc failed!\n");
			continue;
		}

		ldltsl<double,4>(w, d, m, m);

		Vec4d face_dcurv(m[0], m[1], m[2], m[3]);

		// Push it back out to each vertex
		for (int j = 0; j < 3; j++) 
		{
			int vj = vi[j];

			Vec4d this_vert_dcurv;

			proj_dcurv(t, b, face_dcurv, pdir1[vj], pdir2[vj], this_vert_dcurv);

			double wt = cornerareas[i][j] / pointareas[vj];

			dcurv[vj] += wt * this_vert_dcurv;
		}
	}

	Surface_mesh::Vertex_property<Vec4d> my_dcurv = src_mesh->vertex_property<Vec4d>("v:dcurv");
	Surface_mesh::Vertex_iterator vit, vend = src_mesh->vertices_end();

	for (vit = src_mesh->vertices_begin(); vit != vend; ++vit){
		uint i = ((Vertex) vit).idx();
		my_dcurv[vit] = dcurv[i];
	}
}

/* Compute the area "belonging" to each vertex or each corner of a triangle 
(defined as Voronoi area restricted to the 1-ring of a vertex, or to the triangle).*/
void Curvature::computePointAreas(Surface_mesh * src_mesh)
{
	// Get from the Surface_mesh everything you need
	Surface_mesh::Vertex_property<Point> points = src_mesh->vertex_property<Point>("v:point");
	Surface_mesh::Face_iterator fit, fend = src_mesh->faces_end();
	Surface_mesh::Vertex_around_face_circulator fvit, fvend;

	cornerareas = std::vector<Vec3d>(src_mesh->n_faces());
	pointareas = std::vector<double>(src_mesh->n_vertices());

	//#pragma omp parallel for
	for (fit = src_mesh->faces_begin(); fit != fend; ++fit) 
	{
		uint vi[3];
		fvit = src_mesh->vertices(fit);

		vi[0] = ((Vertex)fvit).idx(); Point v0 = points[fvit]; ++fvit;
		vi[1]  = ((Vertex)fvit).idx(); Point v1 = points[fvit]; ++fvit;
		vi[2]  = ((Vertex)fvit).idx(); Point v2 = points[fvit]; ++fvit;

		// Edges
		Point e[3] = {v2 - v1, v0 - v2, v1 - v0};

		// Compute corner weights
		double area = 0.5f * cross(e[0], e[1]).norm();

		double squaredNorm_e0 = e[0].norm() * e[0].norm();
		double squaredNorm_e1 = e[1].norm() * e[1].norm();
		double squaredNorm_e2 = e[2].norm() * e[2].norm();

		double l2[3] = { squaredNorm_e0, squaredNorm_e1, squaredNorm_e2};

		double ew[3] = { l2[0] * (l2[1] + l2[2] - l2[0]), l2[1] * (l2[2] + l2[0] - l2[1]), l2[2] * (l2[0] + l2[1] - l2[2]) };

		uint i = ((Surface_mesh::Face) fit).idx();

		if (ew[0] <= 0.0f) {
			cornerareas[i][1] = -0.25f * l2[2] * area /	dot(e[0] , e[2]);
			cornerareas[i][2] = -0.25f * l2[1] * area /	dot(e[0] , e[1]);
			cornerareas[i][0] = area - cornerareas[i][1] -	cornerareas[i][2];
		} else if (ew[1] <= 0.0f) {
			cornerareas[i][2] = -0.25f * l2[0] * area /	dot(e[1] , e[0]);
			cornerareas[i][0] = -0.25f * l2[2] * area /	dot(e[1] , e[2]);
			cornerareas[i][1] = area - cornerareas[i][2] -	cornerareas[i][0];
		} else if (ew[2] <= 0.0f) {
			cornerareas[i][0] = -0.25f * l2[1] * area /	dot(e[2] , e[1]);
			cornerareas[i][1] = -0.25f * l2[0] * area /	dot(e[2] , e[0]);
			cornerareas[i][2] = area - cornerareas[i][0] -	cornerareas[i][1];
		} else {
			double ewscale = 0.5f * area / (ew[0] + ew[1] + ew[2]);
			for (int j = 0; j < 3; j++)
			{
				cornerareas[i][j] = ewscale * (ew[(j+1)%3] + ew[(j+2)%3]);
			}
		}

		//#pragma omp atomic
		pointareas[vi[0]] += cornerareas[i][0];

		//#pragma omp atomic
		pointareas[vi[1]] += cornerareas[i][1];

		//#pragma omp atomic
		pointareas[vi[2]] += cornerareas[i][2];
	}
}
