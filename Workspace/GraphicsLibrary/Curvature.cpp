#include "Curvature.h"

double Curvature::AREA( double pp1[3], double pp2[3], double pp3[3] )
{
	Vec p1(pp1);
	Vec p2(pp2);
	Vec p3(pp3);

	Vec v1 = p2 - p1;
	Vec v2 = p3 - p1;
	Vec n = v1 ^ v2;
	return 0.5 * n.norm();
}

void Curvature::computePrincipalCurvatures2( Mesh * src_mesh )
{
	this->mesh = src_mesh;

	printf("Computing curvature (1st method).."); CreateTimer(timer);

	int vertex_N = mesh->numberOfVertices();

	k_max = Vector<double>(vertex_N, 0);
	k_min = Vector<double>(vertex_N, 0);

	pdir2 = Vector<Vec>(vertex_N);
	pdir1 = Vector<Vec>(vertex_N);

	Vector<Point3D> vertex = mesh->getVertices();
	Vector<Normal> normal = mesh->getNormals();

	mesh->getUmbrellas();

	for(int i = 0; i<vertex_N; i++)
	{
		int nei_N;

		VertexDetail * vd = mesh->vd(i);
		Umbrella u = mesh->getUmbrella(i);

		// Get number of neighbors + neighbors indices
		nei_N = vd->valence();
		Vector<int> nei ;

		for(HalfEdgeSet::iterator h = u.halfEdge.begin(); h != u.halfEdge.end(); h++)
		{
			int j = h->vertex(0);
			if(j == i)	j = h->vertex(1);

			nei.push_back(j);
		}

		if(nei_N == 0 || nei_N < 3){
			continue;
		}

		//generate weights 
		double total_w = 0;
		double* w = new double[nei_N];

		if(mesh->vd(i)->checkIsBorder())
		{
			w[0] = AREA(vertex[i], vertex[nei[0]], vertex[nei[1]]);
			total_w += w[0];

			for(int j = 1; j<nei_N-1; j++)
			{
				w[j] = AREA(vertex[i], vertex[nei[j]], vertex[nei[j-1]])
					+ AREA(vertex[i], vertex[nei[j]], vertex[nei[(j+1)%nei_N]]);
				total_w += w[j];
			}

			w[nei_N-1] = AREA(vertex[i], vertex[nei[nei_N-2]], vertex[nei[nei_N-1]]);
			total_w += w[nei_N-1]; 
		}
		else
		{
			for(int j=0; j<nei_N; j++)
			{
				w[j] = AREA(vertex[i], vertex[nei[j]], vertex[nei[(j+nei_N-1)%nei_N]])
					+ AREA(vertex[i], vertex[nei[j]], vertex[nei[(j+1)%nei_N]]);
				total_w += w[j];
			}
		}

		if(total_w == 0){
			k_max[i] = k_min[i] = 0;
			pdir2[i][0] = pdir2[i][1] = pdir2[i][2] = 0;
			pdir1[i][0] = pdir1[i][1] = pdir1[i][2] = 0;
			continue;
		}

		double M[3][3];
		M[0][0] = M[0][1] = M[0][2] = M[1][0] = M[1][1] = M[1][2] = M[2][0] = M[2][1] = M[2][2] = 0;

		for(int m=0; m<nei_N; m++)
		{
			int j = nei[m];
			Vec v = vertex[j] - vertex[i];

			//directional curvature

			double len2 = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
			if((double)len2 == 0)
				continue;
			double ki = 2.0*(v * normal[i])/len2;

			//tangent 
			Vec t;
			t[0] = (1.0-normal[i][0]*normal[i][0])*v[0] 
			- normal[i][0]*normal[i][1]*v[1] - normal[i][0]*normal[i][2]*v[2];
			t[1] = -normal[i][0]*normal[i][1]*v[0] + (1.0-normal[i][1]*normal[i][1])*v[1] 
			- normal[i][1]*normal[i][2]*v[2];
			t[2] = -normal[i][0]*normal[i][2]*v[0] - normal[i][1]*normal[i][2]*v[1] 
			+ (1.0-normal[i][2]*normal[i][2])*v[2];


			//another approximation of directinal curvature
			//if concave then -1 times
			//if(DOT(normal[j],t) < 0)
			//ki = -ki;

			double len = t.norm();
			if((double)len != 0){
				t[0] /= len;
				t[1] /= len;
				t[2] /= len;

				for(int row=0; row<3; row++)
					for(int column=0; column<3; column++)
						M[row][column] += w[m]*ki*t[row]*t[column];
			}			
		}
		delete[] w;

		//normalize weight
		if(total_w == 0){
			k_max[i] = k_min[i] = 0;
			pdir2[i][0] = pdir2[i][1] = pdir2[i][2] = 0;
			pdir1[i][0] = pdir1[i][1] = pdir1[i][2] = 0;
			continue;
		}
		for(int j=0; j<3; j++)
			for(int k=0; k<3; k++)
				M[j][k] /= total_w;

		double W[3];

		if(normal[i][0] < 0){
			W[0] = (1.0 - normal[i][0])/sqrt(2.0-2.0*normal[i][0]);
			W[1] = -normal[i][1]/sqrt(2.0-2.0*normal[i][0]);
			W[2] = -normal[i][2]/sqrt(2.0-2.0*normal[i][0]);
		}
		else{
			W[0] = (1.0 + normal[i][0])/sqrt(2.0+2.0*normal[i][0]);
			W[1] = normal[i][1]/sqrt(2.0+2.0*normal[i][0]);
			W[2] = normal[i][2]/sqrt(2.0+2.0*normal[i][0]);
		}

		double Q[3][3];

		for(int row=0; row<3; row++)
			for(int column=0; column<3; column++)
				if(row != column)
					Q[row][column] = -2.0*W[row]*W[column];
				else
					Q[row][column] = 1.0 - 2.0*W[row]*W[column];
		double T1[3] = {Q[1][0], Q[1][1], Q[1][2]};
		double T2[3] = {Q[2][0], Q[2][1], Q[2][2]};

		double M2[2][2];
		M2[0][0] = M[0][0]*Q[0][1]*Q[0][1] + 2.0f*M[0][1]*Q[0][1]*Q[1][1]
		+ M[1][1]*Q[1][1]*Q[1][1] + 2.0f*M[0][2]*Q[0][1]*Q[1][2]
		+ 2.0f*M[1][2]*Q[1][1]*Q[1][2] + M[2][2]*Q[1][2]*Q[1][2];

		M2[0][1] = Q[0][2]*(M[0][0]*Q[0][1] + M[0][1]*Q[1][1] + M[0][2]*Q[1][2])
			+ Q[1][2]*(M[0][1]*Q[0][1] + M[1][1]*Q[1][1] + M[1][2]*Q[1][2])
			+ Q[2][2]*(M[0][2]*Q[0][1] + M[1][2]*Q[1][1] + M[2][2]*Q[1][2]);
		M2[1][0] = M2[0][1];

		M2[1][1] = M[0][0]*Q[0][2]*Q[0][2] + 2.0f*M[0][1]*Q[0][2]*Q[1][2] 
		+ M[1][1]*Q[1][2]*Q[1][2] + 2.0f*M[0][2]*Q[0][2]*Q[2][2]
		+ M[2][2]*Q[2][2]*Q[2][2] + 2.0f*M[1][2]*Q[1][2]*Q[2][2];

		double th = 0;
		if(M2[0][1] != 0)
			th = 0.5*(M2[1][1]-M2[0][0])/M2[0][1];
		double t = 1.0/(fabs(th)+sqrt(th*th+1));
		if(th < 0)
			t = -t;
		double cos = 1.0/sqrt(t*t+1);
		double sin = t*cos;

		double K1_tmp = M2[0][0] - t*M2[0][1];
		double K2_tmp = M2[1][1] + t*M2[0][1];

		double K1 = 3.0*K1_tmp - K2_tmp;
		double K2 = 3.0*K2_tmp - K1_tmp;

		if(K1 > K2){
			k_max[i] = K1;
			pdir2[i][0] = cos*T1[0] - sin*T2[0];
			pdir2[i][1] = cos*T1[1] - sin*T2[1];
			pdir2[i][2] = cos*T1[2] - sin*T2[2];

			k_min[i] = K2;
			pdir1[i][0] = sin*T1[0] + cos*T2[0];
			pdir1[i][1] = sin*T1[1] + cos*T2[1];
			pdir1[i][2] = sin*T1[2] + cos*T2[2];
		}
		else{
			k_min[i] = K1;
			pdir1[i][0] = cos*T1[0] - sin*T2[0];
			pdir1[i][1] = cos*T1[1] - sin*T2[1];
			pdir1[i][2] = cos*T1[2] - sin*T2[2];

			k_max[i] = K2;
			pdir2[i][0] = sin*T1[0] + cos*T2[0];
			pdir2[i][1] = sin*T1[1] + cos*T2[1];
			pdir2[i][2] = sin*T1[2] + cos*T2[2];
		}
	}

	for(int i=0; i<vertex_N; i++)
	{
		Vec n = normal[i] ^ pdir2[i];

		if((n * pdir1[i]) < 0){
			pdir1[i][0] = -pdir1[i][0];
			pdir1[i][1] = -pdir1[i][1];
			pdir1[i][2] = -pdir1[i][2];
		}	
	}

	printf("Done (%d ms).\n", (int)timer.elapsed());
}

// Rotate a coordinate system to be perpendicular to the given normal
void Curvature::rot_coord_sys( const Vec &old_u, const Vec &old_v, const Vec &new_norm, Vec &new_u, Vec &new_v )
{
	new_u = old_u;
	new_v = old_v;
	Vec old_norm = old_u ^ old_v;
	double ndot = old_norm * new_norm;
	if (unlikely(ndot <= -1.0f)) {
		new_u = -new_u;
		new_v = -new_v;
		return;
	}
	Vec perp_old = new_norm - ndot * old_norm;
	Vec dperp = 1.0f / (1 + ndot) * (old_norm + new_norm);
	new_u -= dperp * (new_u * perp_old);
	new_v -= dperp * (new_v * perp_old);
}

// Re-project a curvature tensor from the basis spanned by old_u and old_v
// (which are assumed to be unit-length and perpendicular) to the
// new_u, new_v basis.
void Curvature::proj_curv( const Vec &old_u, const Vec &old_v, double old_ku, double old_kuv, double old_kv, 
						  const Vec &new_u, const Vec &new_v, double &new_ku, double &new_kuv, double &new_kv )
{
	Vec r_new_u, r_new_v;
	rot_coord_sys(new_u, new_v, old_u ^ old_v, r_new_u, r_new_v);

	double u1 = r_new_u * old_u;
	double v1 = r_new_u * old_v;
	double u2 = r_new_v * old_u;
	double v2 = r_new_v * old_v;

	new_ku = old_ku * u1*u1 + old_kuv * (2.0f * u1*v1) + old_kv * v1*v1;
	new_kuv = old_ku * u1*u2 + old_kuv * (u1*v2 + u2*v1) + old_kv * v1*v2;
	new_kv = old_ku * u2*u2 + old_kuv * (2.0f * u2*v2) + old_kv * v2*v2;
}

// And for derivatives of curvature
void Curvature::proj_dcurv(const Vec &old_u, const Vec &old_v, const Vec4 old_dcurv, 
						   const Vec &new_u, const Vec &new_v, Vec4 & new_dcurv)
{
	Vec r_new_u, r_new_v;
	rot_coord_sys(new_u, new_v, old_u ^ old_v, r_new_u, r_new_v);

	double u1 = r_new_u * old_u;
	double v1 = r_new_u * old_v;
	double u2 = r_new_v * old_u;
	double v2 = r_new_v * old_v;

	new_dcurv[0] = old_dcurv[0]*u1*u1*u1 + old_dcurv[1]*3.0f*u1*u1*v1 + old_dcurv[2]*3.0f*u1*v1*v1 + old_dcurv[3]*v1*v1*v1;
	new_dcurv[1] = old_dcurv[0]*u1*u1*u2 + old_dcurv[1]*(u1*u1*v2 + 2.0f*u2*u1*v1) + old_dcurv[2]*(u2*v1*v1 + 2.0f*u1*v1*v2) + old_dcurv[3]*v1*v1*v2;
	new_dcurv[2] = old_dcurv[0]*u1*u2*u2 + old_dcurv[1]*(u2*u2*v1 + 2.0f*u1*u2*v2) + old_dcurv[2]*(u1*v2*v2 + 2.0f*u2*v2*v1) + old_dcurv[3]*v1*v2*v2;
	new_dcurv[3] = old_dcurv[0]*u2*u2*u2 + old_dcurv[1]*3.0f*u2*u2*v2 + old_dcurv[2]*3.0f*u2*v2*v2 + old_dcurv[3]*v2*v2*v2;
}

// Given a curvature tensor, find principal directions and curvatures
// Makes sure that pdir1 and pdir2 are perpendicular to normal
void Curvature::diagonalize_curv( const Vec &old_u, const Vec &old_v, double ku, double kuv, double kv, 
								 const Vec &new_norm, Vec &pdir1, Vec &pdir2, double &k1, double &k2 )
{
	Vec r_old_u, r_old_v;
	rot_coord_sys(old_u, old_v, new_norm, r_old_u, r_old_v);

	double c = 1, s = 0, tt = 0;

	if (likely(kuv != 0.0f)) 
	{
		// Jacobi rotation to diagonalize
		double h = 0.5f * (kv - ku) / kuv;
		tt = (h < 0.0f) ?
			1.0f / (h - sqrt(1.0f + h*h)) :
		1.0f / (h + sqrt(1.0f + h*h));
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

	pdir2 = new_norm ^ pdir1;
}

void Curvature::computePrincipalCurvatures( Mesh * src_mesh )
{
	printf("Computing curvatures (2nd method).."); CreateTimer(timer);

	computePointAreas(src_mesh);

	// Resize the arrays we'll be using
	int nv = vertices.size(), nf = faces.size();
	curv1.clear(); curv1.resize(nv); curv2.clear(); curv2.resize(nv);
	pdir1.clear(); pdir1.resize(nv); pdir2.clear(); pdir2.resize(nv);
	vector<double> curv12(nv);

	// Set up an initial coordinate system per vertex
	for (int i = 0; i < nf; i++) 
	{
		pdir1[faces[i]->VIndex(0)] = vertices[faces[i]->VIndex(1)] - vertices[faces[i]->VIndex(0)];
		pdir1[faces[i]->VIndex(1)] = vertices[faces[i]->VIndex(2)] - vertices[faces[i]->VIndex(1)];
		pdir1[faces[i]->VIndex(2)] = vertices[faces[i]->VIndex(0)] - vertices[faces[i]->VIndex(2)];
	}

	#pragma omp parallel for
	for (int i = 0; i < nv; i++) {
		pdir1[i] = pdir1[i] ^ normals[i];
		pdir1[i].normalize();
		pdir2[i] = normals[i] ^ pdir1[i];
	}

	// Compute curvature per-face
	#pragma omp parallel for
	for (int i = 0; i < nf; i++) {
		// Edges
		Vec e[3] = { vertices[faces[i]->VIndex(2)] - vertices[faces[i]->VIndex(1)],
			vertices[faces[i]->VIndex(0)] - vertices[faces[i]->VIndex(2)],
			vertices[faces[i]->VIndex(1)] - vertices[faces[i]->VIndex(0)] };

		// N-T-B coordinate system per face
		Vec t = e[0];
		t.normalize();
		Vec n = e[0] ^ e[1];
		Vec b = n ^ t;
		b.normalize();

		// Estimate curvature based on variation of normals
		// along edges
		double m[3] = { 0, 0, 0 };
		double w[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };
		for (int j = 0; j < 3; j++) {
			double u = e[j] * t;
			double v = e[j] * b;
			w[0][0] += u*u;
			w[0][1] += u*v;
			//w[1][1] += v*v + u*u; 
			//w[1][2] += u*v; 
			w[2][2] += v*v;
			Vec dn = normals[faces[i]->VIndex(PREV_Index(j))] - normals[faces[i]->VIndex(NEXT_Index(j))];
			double dnu = dn * t;
			double dnv = dn * b;
			m[0] += dnu*u;
			m[1] += dnu*v + dnv*u;
			m[2] += dnv*v;
		}
		w[1][1] = w[0][0] + w[2][2];
		w[1][2] = w[0][1];

		// Least squares solution
		double diag[3];
		if (!ldltdc<double,3>(w, diag)) {
			//dprintf("ldltdc failed!\n");
			continue;
		}
		ldltsl<double,3>(w, diag, m, m);

		// Push it back out to the vertices
		for (int j = 0; j < 3; j++) 
		{
			int vj = faces[i]->VIndex(j);
			double c1, c12, c2;
			proj_curv(t, b, m[0], m[1], m[2], pdir1[vj], pdir2[vj], c1, c12, c2);
			double wt = cornerareas[i][j] / pointareas[vj];

			#pragma omp atomic
			curv1[vj] += wt * c1;

			#pragma omp atomic
			curv12[vj] += wt * c12;

			#pragma omp atomic
			curv2[vj] += wt * c2;
		}
	}

	// Compute principal directions and curvatures at each vertex
	#pragma omp parallel for
	for (int i = 0; i < nv; i++) 
	{
		diagonalize_curv(pdir1[i], pdir2[i],curv1[i], curv12[i], curv2[i], normals[i], pdir1[i], pdir2[i],curv1[i], curv2[i]);
	}

	printf("Done (%d ms).\n", (int)timer.elapsed());
}

// Compute derivatives of curvature
void Curvature::computeDerivativesCurvatures(Mesh * src_mesh)
{
	// Compute principal curvatures and directions
	computePrincipalCurvatures(src_mesh);

	printf("Computing dcurv... ");

	// Resize the arrays we'll be using
	int nv = vertices.size(), nf = faces.size();
	dcurv.clear(); dcurv.resize(nv);

	// Compute derivatives of curvature per-face
	#pragma omp parallel for
	for (int i = 0; i < nf; i++) 
	{
		// Edges
		Vec e[3] = { vertices[faces[i]->VIndex(2)] - vertices[faces[i]->VIndex(1)],
					vertices[faces[i]->VIndex(0)] - vertices[faces[i]->VIndex(2)],
					vertices[faces[i]->VIndex(1)] - vertices[faces[i]->VIndex(0)] };

		// N-T-B coordinate system per face
		Vec t = e[0];
		t.normalize();
		Vec n = e[0] ^ e[1];
		Vec b = n ^ t;
		b.normalize();

		// Project curvature tensor from each vertex into this
		// face's coordinate system
		Vec fcurv[3];
		for (int j = 0; j < 3; j++) 
		{
			int vj = faces[i]->VIndex(j);
			proj_curv(pdir1[vj], pdir2[vj], curv1[vj], 0, curv2[vj],t, b, fcurv[j][0], fcurv[j][1], fcurv[j][2]);

		}

		// Estimate derivatives of curvature based on variation of curvature along edges
		double m[4] = { 0, 0, 0, 0 };
		double w[4][4] = { {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0} };

		for (int j = 0; j < 3; j++) 
		{
			// Variation of curvature along each edge
			Vec dfcurv = fcurv[PREV_Index(j)] - fcurv[NEXT_Index(j)];
			double u = e[j] * t;
			double v = e[j] * b;
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
			//dprintf("ldltdc failed!\n");
			continue;
		}

		ldltsl<double,4>(w, d, m, m);

		Vec4 face_dcurv(m);

		// Push it back out to each vertex
		for (int j = 0; j < 3; j++) 
		{
			int vj = faces[i]->VIndex(j);

			Vec4 this_vert_dcurv;

			proj_dcurv(t, b, face_dcurv, pdir1[vj], pdir2[vj], this_vert_dcurv);

			double wt = cornerareas[i][j] / pointareas[vj];

			dcurv[vj] += wt * this_vert_dcurv;
		}
	}

	printf("Done.\n");
}

/* Compute the area "belonging" to each vertex or each corner
of a triangle (defined as Voronoi area restricted to the 1-ring of
a vertex, or to the triangle).*/
void Curvature::computePointAreas(Mesh * src_mesh)
{
	this->mesh = src_mesh;

	// Get from the mesh everything you need
	vertices = mesh->vertex;
	normals = mesh->vNormal;
	faces = mesh->faceIndexMap;

	printf("Computing point areas... ");

	int nf = faces.size(), nv = vertices.size();
	pointareas.clear();
	pointareas.resize(nv);
	cornerareas.clear();
	cornerareas.resize(nf);

	#pragma omp parallel for
	for (int i = 0; i < nf; i++) 
	{
		// Edges
		Vec e[3] = { vertices[faces[i]->VIndex(2)] - vertices[faces[i]->VIndex(1)],
			vertices[faces[i]->VIndex(0)] - vertices[faces[i]->VIndex(2)],
			vertices[faces[i]->VIndex(1)] - vertices[faces[i]->VIndex(0)] };

		// Compute corner weights
		double area = 0.5f * (e[0] ^ e[1]).norm();
		double l2[3] = { e[0].squaredNorm(), e[1].squaredNorm(), e[2].squaredNorm() };

		double ew[3] = { l2[0] * (l2[1] + l2[2] - l2[0]),
			l2[1] * (l2[2] + l2[0] - l2[1]),
			l2[2] * (l2[0] + l2[1] - l2[2]) };
		if (ew[0] <= 0.0f) {
			cornerareas[i][1] = -0.25f * l2[2] * area /
				(e[0] * e[2]);
			cornerareas[i][2] = -0.25f * l2[1] * area /
				(e[0] * e[1]);
			cornerareas[i][0] = area - cornerareas[i][1] -
				cornerareas[i][2];
		} else if (ew[1] <= 0.0f) {
			cornerareas[i][2] = -0.25f * l2[0] * area /
				(e[1] * e[0]);
			cornerareas[i][0] = -0.25f * l2[2] * area /
				(e[1] * e[2]);
			cornerareas[i][1] = area - cornerareas[i][2] -
				cornerareas[i][0];
		} else if (ew[2] <= 0.0f) {
			cornerareas[i][0] = -0.25f * l2[1] * area /
				(e[2] * e[1]);
			cornerareas[i][1] = -0.25f * l2[0] * area /
				(e[2] * e[0]);
			cornerareas[i][2] = area - cornerareas[i][0] -
				cornerareas[i][1];
		} else {
			double ewscale = 0.5f * area / (ew[0] + ew[1] + ew[2]);
			for (int j = 0; j < 3; j++)
				cornerareas[i][j] = ewscale * (ew[(j+1)%3] +
				ew[(j+2)%3]);
		}

		#pragma omp atomic
		pointareas[faces[i]->VIndex(0)] += cornerareas[i][0];

		#pragma omp atomic
		pointareas[faces[i]->VIndex(1)] += cornerareas[i][1];

		#pragma omp atomic
		pointareas[faces[i]->VIndex(2)] += cornerareas[i][2];
	}

	printf("Done.\n");
}

void Curvature::smoothPrincipalCurvatures(int iterations)
{
	printf("Smoothing principal curvatures... ");

	for(int iter = 0; iter < iterations; iter++)
	{
		int n = (int) pdir1.size();

		Vector<Vec> temp_pdir( n );

		for(int i = 0; i < n; i++)
		{
			Vec newDir;
			StdSet<int> adj = mesh->vertexInfo[i].adjacentVertices();

			foreach(int j, adj)
				newDir += pdir1[j];

			if(newDir.norm() < Epsilon)
				newDir = pdir1[i];

			temp_pdir[i] = newDir.unit();
		}

		for(int i = 0; i < n; i++)
		{
			Vec v = mesh->vec(i);

			Vec planeNormal = pdir1[i] ^ pdir2[i];

			temp_pdir[i].projectOnPlane( planeNormal );

			pdir1[i] = temp_pdir[i];
			//pdir2[i] = pdir1[i].orthogonalVec();
		}
	}

	printf("Done.\n");
}
