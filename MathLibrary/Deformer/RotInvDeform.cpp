#include "RotInvDeform.h"
#include "Utility/SimpleDraw.h"

typedef SparseMatrix<double,RowMajor> SparseMatrixType;

RotInvDeform::RotInvDeform( QSurfaceMesh * usingMesh )
{
	this->mesh = usingMesh;

	// Frequently used
	points = mesh->vertex_property<Point>("v:point");
	normals = mesh->vertex_property<Normal>("v:normal");
	vend = mesh->vertices_end();
	nVerts = mesh->n_vertices();

	globalScale = 1.0;

	ComputeWeights();
	AddBoundaryConstraints();
}

void RotInvDeform::AddBoundaryConstraints( double fWeight /*= 1.0*/ )
{
	for (vit = mesh->vertices_begin(); vit != vend; ++vit)
	{
		if(mesh->is_boundary(vit))
		{
			UpdatePositionConstraint(vit, points[vit], fWeight);
			UpdateOrientationConstraint(vit, vVertices[Surface_mesh::Vertex(vit).idx()].vFrame, fWeight);
		}
	}
}

void RotInvDeform::UpdatePositionConstraint( Surface_mesh::Vertex vID, const Vector3 & vPosition, double fWeight )
{
	bool bFound = false;

	// Update previously set position constraint
	int nCount = vPosConstraints.size();
	for ( int k = 0; !bFound && k < nCount; k++ ) 
	{
		if ( vPosConstraints[k].vID == vID ) 
		{
			vPosConstraints[k].vPosition = vPosition;
			if ( vPosConstraints[k].fWeight != fWeight )
				isMatrixValid = false;
			vPosConstraints[k].fWeight = fWeight;
			bFound = true;
		}
	}

	// Set weight for a new constraint
	if ( ! bFound ) {
		PosConstraint c;
		c.vID = vID;
		c.vPosition = vPosition;
		c.fWeight = fWeight;
		vPosConstraints.push_back(c);
		isMatrixValid = false;
	}
}

void RotInvDeform::UpdateOrientationConstraint( Surface_mesh::Vertex vID, const Frame3 & vFrame, double fWeight )
{
	bool bFound = false;

	// Update previously set orientation constraint
	int nCount = vRotConstraints.size();
	for ( int k = 0; !bFound && k < nCount; k++ ) 
	{
		if ( vRotConstraints[k].vID == vID ){
			vRotConstraints[k].vFrame = vFrame;
			if ( vRotConstraints[k].fWeight != fWeight )
				isMatrixValid = false;
			vRotConstraints[k].fWeight = fWeight;
			bFound = true;
		}
	}

	// Set weight for a new constraint
	if ( ! bFound ) {
		RotConstraint c;
		c.vID = vID;
		c.vFrame = vFrame;
		c.fWeight = fWeight;
		vRotConstraints.push_back(c);
		isMatrixValid = false;
	}
}

Frame3 RotInvDeform::GetCurrentFrame( Surface_mesh::Vertex vID )
{
	Vector3 vNormal = normals[vit];

	// compute frame
	const VtxInfo & vi = vVertices[vID.idx()];
	Vector3 vVtx = points[vID], vNbr = points[vi.nTangentNbr];

	Vector3 xik = vNbr - vVtx;
	Vector3 xikbar = cross(vNormal, cross(xik, vNormal) );
	xikbar.normalize();
	xik = cross( vNormal, xikbar );

	Frame3 vFrame(vVtx);
	vFrame.SetFrame( xikbar, xik, vNormal );
	return vFrame;
}

Frame3 RotInvDeform::GetSolvedFrame( Surface_mesh::Vertex vID )
{
	return vVertices[vID.idx()].vTransFrame;
}

void RotInvDeform::SetVertexWeight( Surface_mesh::Vertex vID, double fWeight )
{
	VtxInfo & vi = vVertices[vID.idx()];
	vi.fVertexWeight = fWeight;
	isMatrixValid = false;
}

void RotInvDeform::SetVertexScale( Surface_mesh::Vertex vID, double fScale )
{
	VtxInfo & vi = vVertices[vID.idx()];
	vi.fVertexScale = fScale;
}

void RotInvDeform::ComputeWeights()
{
	printf("Computing weights..");
	QElapsedTimer timer; timer.start();

	// Clear variables
	vVertices.clear();
	vVertices.resize( mesh->n_vertices() );
	nEdges = 0;
	avgVtxArea = 0.0;

	for (vit = mesh->vertices_begin(); vit != vend; ++vit)
	{
		int vidx = Surface_mesh::Vertex(vit).idx();

		VtxInfo & vi = vVertices[vidx];
		vi.vID = vit;

		// Get one ring of vertex
		vi.vNbrs = mesh->verticesAroundVertex(vit);
		int nNbrs = mesh->valence(vit);
		nEdges += nNbrs;

		// Compute weights for vertex
		CotangentWeights(vit, vi.vNbrWeights);

		Vector3 vVtx, vNormal, vNbr;

		vVtx = points[vit];
		vNormal = normals[vit];

		// find most-orthogonal outgoing edge (for tangent-frame calculation)
		double fMinNbrDot = 1.0;   int nBestNbr = -1;
		for ( int k = 0; k < nNbrs; ++k ) {
			vNbr = points[vi.vNbrs[k]];
			vNbr -= vVtx;   
			vNbr.normalize();
			double fDot = abs(dot(vNbr,vNormal));
			if ( fDot < fMinNbrDot ) {
				fMinNbrDot = fDot;
				nBestNbr = k;
			}
		}
		if ( nBestNbr == -1 ) nBestNbr = 0;
		vi.nTangentNbr = vi.vNbrs[nBestNbr];

		// compute frame
		vNbr = points[vi.nTangentNbr];
		Vector3 xik = vNbr - vVtx;
		Vector3 xikbar = cross(vNormal, cross(xik,vNormal) );
		xikbar.normalize();
		xik = cross(vNormal, xikbar);

		vi.vFrame.SetFrame( xikbar, xik, vNormal);

		//std::cout << "Frame:\n" << vi.vFrame.FrameMatrix();

		// compute laplacian vector
		vi.vLaplacian = MeshLaplacian(vit, vi.vNbrWeights);
		vi.vFrameLaplacian = vi.vLaplacian;
		vi.vFrame.ToFrameLocal( vi.vFrameLaplacian );

		// compute vertex area
		vi.fVertexArea = 1.0;
		//vi.fVertexArea = MeshUtils::VertexArea_Mixed(*Mesh, vit);
		avgVtxArea += (vi.fVertexArea / (double)mesh->n_vertices());

		// additional weighting term
		vi.fVertexWeight = 1.0;

		// per-vertex scaling factor
		vi.fVertexScale = 1.0;
	}

	// rescale vertex areas, so that constraint weight scales are not affected
	for (vit = mesh->vertices_begin(); vit != vend; ++vit)
	{
		VtxInfo & vi = vVertices[Surface_mesh::Vertex(vit).idx()];
		vi.fVertexArea /= avgVtxArea;
	}

	printf("Done. (%d ms).\n", timer.elapsed());
}

void RotInvDeform::CotangentWeights( Surface_mesh::Vertex v, std::vector<double> & vWeights, bool bnormalize /*= false */ )
{
	Vector3 vi = points[v], vj, vo;

	vWeights.resize(mesh->valence(v), 0);
	double fWeightSum = 0;

	// Go around vertex
	int j = 0;
	Surface_mesh::Halfedge_around_vertex_circulator h = mesh->halfedges(v), hend = mesh->halfedges(v);
	do{
		Point vj = points[mesh->to_vertex(h)];

		double a = 0, b = 0;

		Surface_mesh::Halfedge h1 = mesh->next_halfedge(h);
		Surface_mesh::Halfedge h2 = mesh->next_halfedge(mesh->opposite_halfedge(h));

		// Opposing vertices of edge vi-vj
		Point vo_1 = points[mesh->to_vertex(h1)];
		Point vo_2 = points[mesh->to_vertex(h2)];

		if( mesh->is_boundary(mesh->edge(h)) )
		{
			// Give me non-half edge of edge
			Surface_mesh::Halfedge hx = (mesh->is_boundary(mesh->next_halfedge(h)))? 
				mesh->next_halfedge(mesh->opposite_halfedge(h)) : mesh->next_halfedge(h);

			Point vo = points[mesh->to_vertex(hx)];
			a = VectorCot(vi-vo, vj-vo);
		}
		else
		{
			a = VectorCot(vi-vo_1, vj-vo_1);
			b = VectorCot(vi-vo_2, vj-vo_2);
		}

		double dCotSum = a + b;

		vWeights[j++] = dCotSum / 2.0;
		fWeightSum += dCotSum / 2.0;

	} while( ++h != hend );

	if ( bnormalize ) {
		for ( int k = 0; k < mesh->valence(v); ++k )
			vWeights[k] /= fWeightSum;
	}
}

double RotInvDeform::VectorCot( const Vector3 & v1, const Vector3 & v2 )
{
	double fDot = dot(v1,v2);
	return fDot / sqrt( dot(v1, v1) * dot(v2,v2) - fDot*fDot );
}

Vector3 RotInvDeform::MeshLaplacian( Surface_mesh::Vertex v, std::vector<double> & vWeights )
{
	Vector3 vLaplacian(0.0);

	Vector3 vCenter, vNbr;
	vCenter = points[v];

	int i = 0;

	Surface_mesh::Halfedge_around_vertex_circulator h = mesh->halfedges(v), hend = mesh->halfedges(v);
	do{
		Point vj = points[mesh->to_vertex(h)];
		vLaplacian += vWeights[i++] * ( vj - vCenter );
	} while( ++h != hend );

	return vLaplacian;
}

double RotInvDeform::GetLaplacianError()
{
	double fErr = 0.0;

	for ( int i = 0; i < nVerts; ++i )
	{
		VtxInfo & vi = vVertices[i];
		Vector3 vCurLaplacian = MeshLaplacian(vi.vID, vi.vNbrWeights);
		Vector3 vWantLaplacian = vi.vLaplacian * globalScale * vi.fVertexScale;
		//fErr += (vCurLaplacian-vWantLaplacian).SquaredLength();
		vCurLaplacian.normalize();
		vWantLaplacian.normalize();
		double fDot = dot(vCurLaplacian, vWantLaplacian);
		fErr += (1.0 - fDot);
	}
	return fErr;
}


void RotInvDeform::UpdateMatrices()
{
	if ( isMatrixValid )
		return;

	printf("Updating matrices..");
	QElapsedTimer timer; timer.start();

	UpdateRotMatrix();
	UpdatePosMatrix();

	printf("Done. (%d ms).\n", timer.elapsed());

	isMatrixValid = true;
}

void RotInvDeform::UpdateRotMatrix()
{
	//
	//  build system matrix for orientations
	// 

	SparseMatrixType Rs(3*nEdges, 3*nVerts);
	Rs.reserve( (3*nEdges) * 4 );

	int rij = 0;
	for ( int ri = 0; ri < nVerts; ri++ ) 
	{
		VtxInfo & vi = vVertices[ri];
		Surface_mesh::Vertex i = vi.vID;
		int ci = 3*i.idx();
		double fVtxWeight = 1.0 / ( vi.fVertexArea ) / vi.fVertexWeight;

		int nNbrs = vi.vNbrs.size();
		for ( int ni = 0; ni < nNbrs; ++ni ) 
		{
			Surface_mesh::Vertex j = vi.vNbrs[ni];
			int cj = 3*j.idx();
			double fWeight = vi.vNbrWeights[ni];

			Matrix3d Fi = vVertices[i.idx()].vFrame.FrameMatrix();
			Matrix3d Fj = vVertices[j.idx()].vFrame.FrameMatrix();
			Matrix3d Rij = Fj.transpose() * Fi;

			for ( int k = 0; k < 3; ++k ) 
			{
				Rs.insert( rij+k, ci+0 ) = fWeight * Rij(k, 0);
				Rs.insert( rij+k, ci+1 ) = fWeight * Rij(k, 1);
				Rs.insert( rij+k, ci+2 ) = fWeight * Rij(k, 2);

				Rs.insert( rij+k, cj+k ) = fWeight * -1.0;
			}
			rij += 3;
		}
	}
	SparseMatrixType RsT_Rs = SparseMatrixType(Rs.transpose()) * Rs;


	RotA.resize(3 * nVerts, 3 * nVerts);
	for (int k = 0; k < 3*nVerts; ++k) {
		for ( SparseMatrixType::InnerIterator it(RsT_Rs,k); it; ++it)
			RotA.insert( it.row(), it.col() ) = it.value();
	}

	// Rotation constraints
	int nRotCons = vRotConstraints.size();
	for ( unsigned int ci = 0; ci < nRotCons; ++ci ) 
	{
		RotConstraint & c = vRotConstraints[ci];
		int ri = c.vID.idx() * 3;

		for ( int k = 0; k < 3; ++k )
			RotA.coeffRef(ri+k,ri+k) = RotA.coeffRef(ri+k,ri+k) + c.fWeight*c.fWeight;
	}	

	//IOFormat OctaveFmt(4, 0, ", ", ";\n", "", "", "[", "]");
	//std::cout << "Rotation System:\n" << RotA.toDense().format(OctaveFmt);
}

void RotInvDeform::UpdatePosMatrix()
{
	//
	// build Laplacian system to find positions
	//

	SparseMatrixType Ls(nVerts,nVerts);
	for ( int ri = 0; ri < nVerts; ++ri ) {
		VtxInfo & vi = vVertices[ri];
		int nNbrs = vi.vNbrs.size();

		double dSum = 0.0;
		for ( int k = 0; k < nNbrs; ++k ) {
			Ls.insert(ri, vi.vNbrs[k].idx()) = vi.vNbrWeights[k];
			dSum += vi.vNbrWeights[k];
		}
		Ls.insert(ri, ri) = -dSum;
	}
	SparseMatrixType LsLs = Ls * Ls;


	PosA.resize(nVerts,nVerts);
	for (unsigned int k=0; k < nVerts; ++k) {
		for ( SparseMatrixType::InnerIterator it(LsLs,k); it; ++it)
			PosA.insert(it.row(), it.col()) = it.value();
	}

	// Position constraints
	int nCons = vPosConstraints.size();
	for ( unsigned int ci = 0; ci < nCons; ++ci ) {
		PosConstraint & c = vPosConstraints[ci];
		int ri = c.vID.idx();
		PosA.coeffRef(ri,ri) = PosA.coeffRef(ri,ri) + c.fWeight*c.fWeight;
	}

	//IOFormat OctaveFmt(4, 0, ", ", ";\n", "", "", "[", "]");
	//std::cout << "\nPosition system:\n" << PosA.toDense().format(OctaveFmt);
}

void RotInvDeform::UpdateRHSRot()
{
	RHSRot = MatrixXd::Zero(3 * nVerts, 3);

	// update rotation constraints
	int nRotCons = vRotConstraints.size();
	for ( int ci = 0; ci < nRotCons; ++ci ) 
	{
		RotConstraint & c = vRotConstraints[ci];
		int ri = c.vID.idx() * 3;

		Vector3 vConsValA = c.fWeight * c.fWeight*c.vFrame.X();
		RHSRot.row(ri) = V2E(vConsValA);

		Vector3 vConsValB = c.fWeight * c.fWeight*c.vFrame.Y();
		RHSRot.row(ri+1) = V2E(vConsValB);

		Vector3 vConsValN = c.fWeight * c.fWeight*c.vFrame.Z();
		RHSRot.row(ri+2) = V2E(vConsValN);
	}

	//std::cout << "\n RHS \n" << RHSRot;
}

void RotInvDeform::UpdateRHSPos()
{
	// After solving for orientation,

	// transform laplacians to new frames
	for ( int ri = 0; ri < nVerts; ++ri ) 
	{
		VtxInfo & vi = vVertices[ri];

		// transform frame-encoded laplacian vector into new 3D frame
		vi.vLaplacian = vi.vFrameLaplacian;
		vi.vTransFrame.ToWorld( vi.vLaplacian );
	}

	RHSPos.resize(3, nVerts);

	double fGlobalScale = globalScale;
	for ( int ri = 0; ri < nVerts; ri++ ) 
	{
		VtxInfo & vi = vVertices[ri];
		double x[3] = {0,0,0};
		double fWeightSum = 0.0;
		int nNbrs = vi.vNbrs.size();

		for ( int k = 0; k < nNbrs; ++k ) {
			VtxInfo & ni = vVertices[ vi.vNbrs[k].idx() ];
			double fWeight = vi.vNbrWeights[k];
			fWeightSum += fWeight;
			const Vector3 & vNbrLaplacian = ni.vLaplacian;
			for ( int i = 0; i < 3; ++i ) {
				x[i] += fWeight*vNbrLaplacian[i]*fGlobalScale*ni.fVertexScale;
			}
		}

		const Vector3 & vLaplacian = vi.vLaplacian;
		for ( int i = 0; i < 3; ++i )
			RHSPos.coeffRef(i,ri) = x[i] - fWeightSum*vLaplacian[i]*fGlobalScale*vi.fVertexScale ;
	}

	// add position constraints
	int nCons = vPosConstraints.size();
	for ( int ci = 0; ci < nCons; ++ci ) 
	{
		PosConstraint & c = vPosConstraints[ci];
		int ri = c.vID.idx();
		double fConsWeight = c.fWeight;
		for ( int k = 0; k < 3; ++k ) 
			RHSPos.coeffRef(k,ri) += c.vPosition[k]*fConsWeight*fConsWeight;
	}

	// DEBUG:
	/*std::cout << "\nPosRHS\n";
	for ( int ri = 0; ri < nVerts; ri++ ) 
		std::cout << RHSPos[0][ri] << "\t" << RHSPos[1][ri] << "\t" << RHSPos[2][ri] << "\n";
	int x = 0;*/
}


void RotInvDeform::Solve()
{
	UpdateMatrices();

	printf("Solving..");
	QElapsedTimer timer; timer.start();

	UpdateRHSRot();

	SimplicialLLT< SparseMatrix<double> > rot_solver(RotA);
	std::vector<VectorXd> rot(3, VectorXd::Zero(nVerts*3));

	// Solve for each basis
	for(int i = 0; i < 3; i++)
	{
		VectorXd b = RHSRot.col(i);
		rot[i] = rot_solver.solve( b );
	}
	
	//bool bOKRot = rot_solver.info() == ComputationInfo::Success;
	
	// extract solved frames 
	for ( int i = 0; i < nVerts; ++i )
	{
		VtxInfo & vi = vVertices[i];
		int ri = 3*i;
		Vector3 vFrameV[3];
		for ( int k = 0; k < 3; ++k ) {
			vFrameV[k] = Vector3( rot[0][ri+k], rot[1][ri+k], rot[2][ri+k] );
			vFrameV[k].normalize();
		}
		Vector3 vA = cross(vFrameV[1],vFrameV[2]);		vA.normalize();
		Vector3 vB = cross(vFrameV[2], vA );			vB.normalize();
		vi.vTransFrame.SetFrame( vA, vB, vFrameV[2] );
	}

	UpdateRHSPos();
	
	// Solve
	SimplicialLLT< SparseMatrix<double> > pos_solver(PosA);
	std::vector<VectorXd> pos(3, VectorXd::Zero(nVerts));

	for(int i = 0; i < 3; i++)
	{
		VectorXd b = RHSPos.row(i);
		pos[i] = pos_solver.solve( b );
	}

	// extract solved position
	for ( int i = 0; i < nVerts; ++i ) 
	{
		Vector3 v(pos[0][i], pos[1][i], pos[2][i]);
		points[Surface_mesh::Vertex(i)] = v;
	}

	printf("Done. (%d ms).\n", timer.elapsed());
}

void RotInvDeform::draw()
{
	std::vector<Vec3d> cpoint, bpoint;

	foreach(PosConstraint c, vPosConstraints)
	{
		if(mesh->is_boundary(c.vID))
			bpoint.push_back(points[c.vID]);
		else
			cpoint.push_back(points[c.vID]);
	}

	SimpleDraw::IdentifyPoints(cpoint);
	SimpleDraw::IdentifyPoints(bpoint, Vec4d(1,1,0,1));
}
