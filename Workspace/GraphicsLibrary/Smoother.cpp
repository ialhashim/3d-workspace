#include "Smoother.h"


void Smoother::LaplacianSmoothing(Mesh * m, StdSet<Face*> & faceList, int numIteration, bool protectBorders)
{
	CreateTimer(timer);
	printf("\nPerforming Laplacian smoothing (iterations = %d)...", numIteration);

	if(protectBorders)
		m->flagBorderVertices();

	// This method uses the basic equal weights Laplacian operator
	for(int iteration = 0; iteration < numIteration; iteration++)
	{
		Vector<Vec> newPositions( m->numberOfVertices() );
		Vector<int> count( m->numberOfVertices() );

		StdSet<int> modifiedPoints;

		// For each vertex in all faces, add up the neighbors and
		// store in "newPositions[x]", also keep count of neighbors in 
		// "count[x]"
		for(StdSet<Face*>::iterator f = faceList.begin(); f != faceList.end(); f++)
		{
			Face * face = *f;

			for(int j = 0; j < 3; ++j)
			{
				int j2 = (j+1) % 3;

				int v1 = face->VIndex(j);
				int v2 = face->VIndex(j);

				newPositions[v1] += face->vec(j2);
				newPositions[v2] += face->vec(j);

				count[v1]++;
				count[v2]++;

				modifiedPoints.insert(v1);
				modifiedPoints.insert(v2);
			}
		}

		// New vertex position = (x_old + x_new) / (count )
		for(StdSet<int>::iterator it = modifiedPoints.begin(); it != modifiedPoints.end(); it++)
		{
			int vi = *it;

			newPositions[vi] /= (double)(count[vi]);

			if(protectBorders)
			{
				if(m->vertexInfo[vi].flag != VF_BORDER)
					m->vertex[vi].set( newPositions[vi] );
			}
			else
			{
				m->vertex[vi].set( newPositions[vi] );
			}
		}
	}

	printf("done. (%d ms)\n", (int)timer.elapsed());
}

void Smoother::LaplacianSmoothing(Mesh * m, int numIteration, bool protectBorders)
{
	StdSet<Face*> faceList;

	for(StdList<Face>::iterator f = m->face.begin(); f != m->face.end(); f++)
		faceList.insert(&(*f));	

	LaplacianSmoothing(m, faceList, numIteration, protectBorders);

	m->computeNormals();
	m->computeBounds();

	m->vbo->setDirty(true);
}

Point3D Smoother::LaplacianSmoothVertex(Mesh * m, int vi)
{
	Point3D newPos;

	StdSet<int> adj = m->vertexInfo[vi].adjacentVertices();

	foreach(int j, adj)
		newPos += m->vertex[j];

	newPos /= (double)adj.size();

	return newPos;
}

Point3D Smoother::ScaleDependentSmoothVertex(Mesh * m, int vi, float step_size)
{
	Point3D pt = m->vertex[vi];

	Umbrella * u, tempU;

	if(m->tempUmbrellas.size() < 1)
		u = &m->tempUmbrellas[vi];
	else
	{
		tempU =  m->getUmbrella(vi);
		u = &tempU;
	}

	Vector<Vec> nighbour;
	Vector<double> weight;

	Point3D diff, sigma;
	double E = 0.0;

	for(HalfEdgeSet::iterator h = u->halfEdge.begin(); h != u->halfEdge.end(); h++)
	{
		int j = h->vertex(0);

		if(j == vi)	j = h->vertex(1);

		diff = m->vertex[j] - pt;

		double edgeLen = diff.norm();

		if(edgeLen == 0){
			E = 0;	break;
		}

		E += edgeLen;

		weight.push_back(edgeLen);
		nighbour.push_back(m->vertex[j]);
	}

	// Normalize weights
	if(E > 0){
		for(int j = 0; j < (int)nighbour.size(); j++)
		{
			weight[j] /= E;
			sigma += (weight[j] * (nighbour[j] - pt));
		}
	}
	else
		return pt;

	return pt + (sigma * step_size);
}

void Smoother::ScaleDependentSmoothing(Mesh * m, int numIteration, float step_size, bool protectBorders)
{
	printf("\nPerforming Scale Dependent Smoothing (iterations = %d)...", numIteration);

	if(protectBorders)
		m->flagBorderVertices();

	Vector<Umbrella> * U;
	if(m->tempUmbrellas.size() < 1)	m->getUmbrellas();
	U = &m->tempUmbrellas;

	for(int iteration = 0; iteration < numIteration; iteration++)
	{
		Vector<Vec> newPositions( m->numberOfVertices() );

		for(int vi = 0; vi < m->numberOfVertices(); vi++)
		{
			Point3D pt = m->vertex[vi];

			if(m->vertexInfo[vi].flag != VF_BORDER)
			{
				Vector<Vec> nighbour;
				Vector<double> weight;

				Point3D diff, sigma;

				double E = 0.0;

				for(HalfEdgeSet::iterator h = U->at(vi).halfEdge.begin(); h != U->at(vi).halfEdge.end(); h++)
				{
					int j = h->vertex(0);
					if(j == vi)	j = h->vertex(1);

					diff = m->vertex[j] - pt;

					double edgeLen = diff.norm();

					if(edgeLen == 0)
					{
						E = 0;
						break;
					}

					E += edgeLen;

					weight.push_back(edgeLen);
					nighbour.push_back(m->vertex[j]);
				}

				// Normalize weights
				if(E > 0)
				{
					for(int j = 0; j < (int)nighbour.size(); j++)
					{
						weight[j] /= E;
						sigma += (weight[j] * (nighbour[j] - pt));
					}

					newPositions[vi] = pt + (sigma * step_size);
				}
			}
		}

		#pragma omp parallel for
		for(int i = 0; i < (int)newPositions.size(); i++)
		{
			if(m->vertexInfo[i].flag != VF_BORDER)
				m->vertex[i].set( newPositions[i] );
		}
	}

	m->computeNormals();
	m->computeBounds();

	m->vbo->setDirty(true);
}

/* 
* Implementation of "Implicit Fairing of Irregular Meshes using Diffusion and Curvature Flow"
* Based on http://www.stanford.edu/~vkl/code/vdec.html
*
* Instead of (I - lambda_dt K) X_{n+1} = X_n we do:
* (A - lambda_dt Ks) X_{n+1} = A X_n  (multiplie both sides by A)
*
* K_{ij} = K_{ji} (hopefully)
* K_{ij} = cot_{ij} for i neq j where cot_{ij} is the cotangent term for edge ij
* K_{ij} = -sum_j cot_{ij} for i = j
*
*/
void Smoother::MeanCurvatureFlow(Mesh * mesh, double step,int numIteration, bool isVolumePreservation)
{
	if(step == 0.0)	return;

        mesh = mesh;
        numIteration = numIteration;
        isVolumePreservation = isVolumePreservation;

	// WARNGING: this is disabled since SparseLib++ is not included in this project (for simplicity)

	/*CreateTimer(timer);
	printf("\n\nPerforming Mean Curvature Flow smoothing (iterations = %d, step = %f)...", numIteration, step);

	mesh->flagBorderVertices();

	int real_N = mesh->numberOfVertices();

	// Should be replaced by HalfEdges ?
	Vector<Umbrella> * U;

	if(mesh->tempUmbrellas.size() < 1)
		mesh->getUmbrellas();

	U = &mesh->tempUmbrellas;

	// Process borders, from paper's suggestion of virtual center point
	treatBorders(mesh, *U);

	int N = U->size();

	// Initialize B & X vectors

	VECTOR_double b_x(N), b_y(N), b_z(N);
	VECTOR_double X(N), Y(N), Z(N);
	Vector<double> area (N, 0.0);

	double dt = step;
	double alpha, beta, cots;

	double init_volume = mesh->computeVolume();

	for(int k = 0; k < numIteration; k++)
	{
		// Initialize dynamic sparse matrix, used for creation of sparse M
		Eigen::DynamicSparseMatrix<double> M(N, N);
		M.reserve(N * 7);

		// For each point
		for(int i = 0; i < N; i++)
		{
			Vector<Face *> * ifaces = &(U->at(i)).ifaces;

			// for each nighboring edge
			for(int f = 0; f < (int)ifaces->size(); f++)
				area[i] += ifaces->at(f)->area();

			if(area[i] == 0) printf(".zero area."); // should not happen

			M.coeffRef(i,i) = area[i];

			if((U->at(i)).flag != VF_BORDER)
			{
				for(HalfEdgeSet::iterator halfEdge = (U->at(i)).halfEdge.begin(); halfEdge != (U->at(i)).halfEdge.end(); halfEdge++)
				{
					int j = halfEdge->vertex(0);
					if(j == i)	j = halfEdge->vertex(1);

					if(halfEdge->face(0) == NULL || halfEdge->face(1) == NULL)
					{
						cots = 0;
					}
					else
					{
						Edge e = halfEdge->edge;
						alpha = halfEdge->face(0)->angleCotangent(e);
						beta = halfEdge->face(1)->angleCotangent(e);

						cots = (alpha + beta) * 0.25 * dt;

						if(alpha == 0 || beta == 0 || (U->at(j)).flag == VF_BORDER)
							cots = 0;
					}

					if(i > j)
					{
						M.coeffRef(i,j) -= cots;
						M.coeffRef(j,i) -= cots;
					}

					M.coeffRef(i,i) += cots;
				}
			}

			Vertex * vertex;

			if( i < real_N )
				vertex = mesh->v(i);
			else
				vertex =  (U->at(i)).ifaces[0]->PointIndexed(i); // Special border virtual vertices

			X(i) = vertex->x;
			Y(i) = vertex->y;
			Z(i) = vertex->z;

			b_x(i) = area[i] * X(i);
			b_y(i) = area[i] * Y(i);
			b_z(i) = area[i] * Z(i);
		}

		Eigen::SparseMatrix<double> mat(M);
		CompCol_Mat_double A(N, N, mat.nonZeros(), mat._valuePtr(), mat._innerIndexPtr(), mat._outerIndexPtr());

		//Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
		//std::cout << "\n" << mat.toDense().format(CleanFmt) << "\n";

		// Solver tolerance
		double tol = 0.01;
		double tol_x = tol, tol_y = tol, tol_z = tol;
		int maxit = 500, maxit_x = maxit, maxit_y = maxit, maxit_z = maxit;
		int result;

		// Precondition
		DiagPreconditioner_double precond(A);

		CreateTimer(solveTime);
		printf("\nSolving..");

		// Solve
		//		#pragma omp parallel sections
		{
			//			#pragma omp section
			{
				result = Solver::BiCG(A, X, b_x, precond, maxit_x, tol_x);
				printf(" conv X = %s ..",(result)?"false":"true");
			}

			//			#pragma omp section
			{
				result = Solver::BiCG(A, Y, b_y, precond, maxit_y, tol_y);
				printf(" conv Y = %s ..",(result)?"false":"true");
			}

			//			#pragma omp section
			{
				result = Solver::BiCG(A, Z, b_z, precond, maxit_z, tol_z);
				printf(" conv Z = %s ..",(result)?"false":"true");
			}
		}

		printf("\nSolve time = %d ms\n", (int)solveTime.elapsed());

		Vec center = mesh->computeCenter();

		//if(result == 0)
		{
			// Assign returned solution
			for(int i = 0; i < real_N; i++)
				mesh->vertex[i].set(X(i), Y(i), Z(i));
		}

		mesh->translate(-center);

		if(isVolumePreservation)
		{
			// Volume preservation
			double new_volume = mesh->computeVolume();
			mesh->scale(pow(init_volume / new_volume, 1.0 / 3.0));
		}

		mesh->translate(center);

		// Undo addition of virtual hole center points
		untreatBorders(mesh, *U);

	}

	printf("Smoothing done. (%d ms)\n", (int)timer.elapsed());

	mesh->computeNormals();
	mesh->computeBounds();

	mesh->vbo->setDirty(true);*/
}

void Smoother::MeanCurvatureFlowExplicit(Mesh * mesh, double step, int numIteration)
{
	CreateTimer(timer);
	printf("\n\nPerforming Mean Curvature Flow smoothing (Explicit, iterations = %d)...", numIteration);
	Vector<Umbrella> * U;

	if(mesh->tempUmbrellas.size() < 1)
		mesh->getUmbrellas();

	U = &mesh->tempUmbrellas;

	// borders code here v

	int N = U->size();

	for(int iteration = 0; iteration < numIteration; iteration++)
	{
		Vector<Vec> pos( mesh->numberOfVertices() );
		Vector<double> cots( mesh->numberOfVertices() );

		double init_volume = mesh->computeVolume();

		double alpha, beta;

		// For each point
		for(int i = 0; i < N; i++)
		{
			if(U->at(i).flag != VF_BORDER)
			{
				for(HalfEdgeSet::iterator halfEdge = U->at(i).halfEdge.begin(); halfEdge != U->at(i).halfEdge.end(); halfEdge++)
				{
					int j = halfEdge->vertex(0);
					if(j == i)	j = halfEdge->vertex(1);

					double cot = 0;

					if(halfEdge->face(0) == NULL || halfEdge->face(1) == NULL || U->at(j).flag == VF_BORDER)
					{
						cot = 0;
					}
					else
					{
						Edge e = halfEdge->edge;
						alpha = halfEdge->face(0)->angleCotangent(e);
						beta = halfEdge->face(1)->angleCotangent(e);

						cot = alpha + beta;

						if(alpha == 0 || beta == 0)
							cot = 0;				

						cots[i] += cot;
					}

					pos[i] += mesh->vertex[j] * cot;
				}
			}
		}

		for(int i = 0; i < (int)pos.size(); i++)
		{
			Vec center = mesh->vertex[i];

			if(cots[i] > 0)
			{
				pos[i] *= (step / cots[i]);
				pos[i] += ((1.0 - step) * center);
			}
			else
			{
				pos[i] = center;
			}
		}

		for(int i = 0; i < (int)pos.size(); i++)
		{
			mesh->vertex[i] = pos[i];
		}

		Vec center = mesh->computeCenter();

		mesh->translate(-center);
		// Volume preservation
		double new_volume = mesh->computeVolume();
		mesh->scale(pow(init_volume / new_volume, 1.0 / 3.0));
		mesh->translate(center);
	}

	printf(" done. (%d ms)\n", (int)timer.elapsed());

	mesh->computeNormals();
	mesh->computeBounds();

	mesh->vbo->setDirty(true);
}

void Smoother::treatBorders(Mesh * mesh, Vector<Umbrella> & U)
{
	// Debug
	//mesh->markers.clear();
	//mesh->testFaces.clear();

	HoleStructure hole = mesh->getHoles();

	// If there is a border, there is a hole
	if(hole.size())
	{
		// Reusable elements
		int v1_index, v2_index, v3_index;
		Point3D *v2, *v3;

		// For each hole, fill with virtual faces
		for(StdMap<int, Vector<int> >::iterator it = hole.begin(); it != hole.end(); it++)
		{
			Vector<int> * points = &it->second;

			// A non-manifold is not a hole
			if(points->size() > 2)
			{
				Point3D sum, *hole_vertex;

				// Find center point
				for(int i = 0; i < (int)points->size(); i++)
				{
					//mesh->markers[mesh->v(points->at(i))] = Color4(255,0,0);
					sum += *mesh->v(points->at(i));
				}

				hole_vertex = new Point3D(sum / points->size());

				VertexDetail * hold_detail = new VertexDetail(U.size());
				U.push_back(Umbrella(hold_detail));

				// Create faces to connect each hole point
				for(int i = 0; i < (int)points->size(); i++)
				{
					v1_index = hold_detail->index;
					v2_index = points->at(i);
					v3_index = points->at((i+1) % points->size());

					v2 = mesh->v(v2_index);
					v3 = mesh->v(v3_index);

					// Create the face
					Face * f = new Face(v1_index, v2_index, v3_index, hole_vertex, v2, v3, hold_detail->index);
					//mesh->testFaces.push_back(*f); // debug

					// Add face to incident structures of involved vertices
					U[v1_index].ifaces.push_back(f);
					U[v2_index].ifaces.push_back(f);
					U[v3_index].ifaces.push_back(f);

					// Unflag as border
					U[v2_index].flag = VF_CLEAR;
					U[v3_index].flag = VF_CLEAR;
				}

				// Find umbrella around virtual point and its connected vertices
				U[hold_detail->index].loadHalfEdgeSet();

				for(int i = 0; i < (int)points->size(); i++)
					U[points->at(i)].loadHalfEdgeSet();
			}
		}
	}
}

void Smoother::untreatBorders(Mesh * mesh, Vector<Umbrella> & U)
{
	int N = mesh->numberOfVertices();

	int numBorders = U.size() - N;

	if(numBorders > 0)
	{
		for(int i = 0 ; i < numBorders; i++)
		{
			int vi = U.size() - i - 1;

			// Refresh modified umbrellas
			Umbrella * u = &U.at(vi);
			foreach(int j, u->neighbor)
			{
				mesh->tempUmbrellas[j] = Umbrella(mesh->vd(j));
			}
		}

		U.resize(N);
	}
}
