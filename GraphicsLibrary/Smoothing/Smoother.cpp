#include "Smoother.h"

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Core>
#include <Eigen/Sparse>
using namespace Eigen;

// SparseLib++ : Linear Solver
#include "SparseLib/sparse_lib.h"

void Smoother::LaplacianSmoothing(Surface_mesh * m, int numIteration, bool protectBorders)
{
	//CreateTimer(timer);
	printf("\nPerforming Laplacian smoothing (iterations = %d)...", numIteration);

	Surface_mesh::Vertex_property<Point> points = m->vertex_property<Point>("v:point");
	Surface_mesh::Vertex_property<Point> newPositions = m->vertex_property<Point>("v:new_point");

	Surface_mesh::Vertex_around_vertex_circulator vvit, vvend;

	// This method uses the basic equal weights Laplacian operator
	for(int iteration = 0; iteration < numIteration; iteration++)
	{
		Surface_mesh::Vertex_iterator vit, vend = m->vertices_end();

		// Original positions, for boundary
		for(vit = m->vertices_begin(); vit != vend; ++vit)
			newPositions[vit] = points[vit];
		
		// Compute Laplacian
		for(vit = m->vertices_begin(); vit != vend; ++vit)
		{
			if(!protectBorders || (protectBorders && !m->is_boundary(vit)))
			{
				newPositions[vit] = Point(0,0,0);

				// Sum up neighbors
				vvit = vvend = m->vertices(vit);
				do{ newPositions[vit] += points[vvit]; } while(++vvit != vvend);

				// Average it
				newPositions[vit] /= m->valence(vit);
			}
		}

		// Set vertices to final position
		for(vit = m->vertices_begin(); vit != vend; ++vit)
			points[vit] = newPositions[vit];
	}

	//printf("done. (%d ms)\n", (int)timer.elapsed());
}

Point Smoother::LaplacianSmoothVertex(Surface_mesh * m, int vi)
{
	Point newPos(0,0,0);

	Surface_mesh::Vertex_property<Point> points = m->vertex_property<Point>("v:point");
	Surface_mesh::Vertex_around_vertex_circulator vit, vend;

	vit = vend = m->vertices(Surface_mesh::Vertex(vi));
	do{ newPos += points[vit]; } while(++vit != vend);

	return newPos;
}

Point Smoother::ScaleDependentSmoothVertex(Surface_mesh * m, Surface_mesh::Vertex v, double step_size)
{
	Surface_mesh::Vertex_property<Point> points = m->vertex_property<Point>("v:point");

	Point pt = points[v];

	std::vector<Point> nighbour;
	std::vector<double> weight;

	Point diff(0,0,0), sigma(0,0,0);
	double E = 0.0;

	Surface_mesh::Halfedge_around_vertex_circulator hit = m->halfedges(v), hend = m->halfedges(v);
	do{
		Surface_mesh::Vertex vj = m->to_vertex(hit);

		double edgeLen = (points[vj] - pt).norm();

		if(edgeLen == 0){
			E = 0;	
			break;
		}

		E += edgeLen;

		weight.push_back(edgeLen);
		nighbour.push_back(points[vj]);

	} while(++hit != hend);

	// Normalize weights
	if(E > 0)
	{
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

void Smoother::ScaleDependentSmoothing(Surface_mesh * m, int numIteration, double step_size, bool protectBorders)
{
	printf("\nPerforming Scale Dependent Smoothing (iterations = %d)...", numIteration);

	Surface_mesh::Vertex_property<Point> points = m->vertex_property<Point>("v:point");
	Surface_mesh::Vertex_property<Point> newPositions = m->vertex_property<Point>("v:new_point");
	Surface_mesh::Vertex_iterator vit, vend = m->vertices_end();

	for(int iteration = 0; iteration < numIteration; iteration++)
	{
		// Original positions, for boundary
		for(vit = m->vertices_begin(); vit != vend; ++vit)
			newPositions[vit] = points[vit];

		// Compute scale-dependent Laplacian
		for(vit = m->vertices_begin(); vit != vend; ++vit)
		{
			if(!protectBorders || (protectBorders && !m->is_boundary(vit)))
				newPositions[vit] = ScaleDependentSmoothVertex(m, vit, step_size);
		}

		// Set vertices to final position
		for(vit = m->vertices_begin(); vit != vend; ++vit)
			points[vit] = newPositions[vit];
	}
}

/* 
* Implementation of "Implicit Fairing of Irregular QSurfaceMeshes using Diffusion and Curvature Flow"
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
void Smoother::MeanCurvatureFlow(QSurfaceMesh * m, int numIteration, double step, bool isVolumePreservation)
{
	if(step == 0.0)	return;

	Surface_mesh::Vertex_property<Point> points = m->vertex_property<Point>("v:point");
	Surface_mesh::Face_around_vertex_circulator fvit, fvend;
	Surface_mesh::Halfedge_around_vertex_circulator h, hend;

	numIteration = numIteration;
	isVolumePreservation = isVolumePreservation;

	// WARNGING: this is disabled since SparseLib++ is not included in this project (for simplicity)

	//CreateTimer(timer);
	printf("\n\nPerforming Mean Curvature Flow smoothing (iterations = %d, step = %f)...", numIteration, step);

	// Process borders, from paper's suggestion of virtual center point
	//treatBorders(QSurfaceMesh, *U);

	int N = m->n_vertices();

	// Initialize B & X vectors
	VECTOR_double b_x(N), b_y(N), b_z(N);
	VECTOR_double X(N), Y(N), Z(N);
	std::vector<double> area (N, 0.0);

	double dt = step;
	double alpha, beta, cots, cot_alpha, cot_beta;

	double init_volume = m->volume();

	// Perform smoothing iteration
	for(int k = 0; k < numIteration; k++)
	{
		// Initialize dynamic sparse matrix, used for creation of sparse M
		Eigen::DynamicSparseMatrix<double> M(N, N);
		M.reserve(N * 7);

		// For each point
		for(int i = 0; i < N; i++)
		{
			Surface_mesh::Vertex v(i);

			// for each neighboring face
			fvit = fvend = m->faces(v);
			do{	area[i] += m->faceArea(fvit); } while (++fvit != fvend);

			M.coeffRef(i,i) = area[i];

			// for each neighboring vertex
			if(!m->is_boundary(v))
			{
				Point p0 = points[v];

				h = m->halfedges(v), hend = m->halfedges(v);
				do{
					Surface_mesh::Vertex vj = m->to_vertex(h);
					int j = vj.idx();

					Point p1 = points[vj];
					Point p2 = points[m->to_vertex(m->next_halfedge(h))];
					Point p3 = points[m->to_vertex(m->next_halfedge(m->opposite_halfedge(h)))];

					double cos_alpha = dot((p0 - p2).normalize(), (p1 - p2).normalize());
					double cos_beta = dot((p0 - p3).normalize(), (p1 - p3).normalize());

					cot_alpha = cos_alpha / sin(acos( cos_alpha ));
					cot_beta = cos_beta / sin(acos( cos_beta ));

					cots = (cot_alpha + cot_beta) * 0.25 * dt;

					if(alpha == 0 || beta == 0)
						cots = 0;
					
					if(i > j)
					{
						M.coeffRef(i,j) -= cots;
						M.coeffRef(j,i) -= cots;
					}

					M.coeffRef(i,i) += cots;

				} while( ++h != hend );
			}

			X(i) = points[v].x();
			Y(i) = points[v].y();
			Z(i) = points[v].z();

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

		//CreateTimer(solveTime);
		printf("\nSolving..");

		// Solve
		#pragma omp parallel sections
		{
			#pragma omp section
			{
				result = Solver::BiCG(A, X, b_x, precond, maxit_x, tol_x);
				printf(" conv X = %s ..",(result)?"false":"true");
			}

			#pragma omp section
			{
				result = Solver::BiCG(A, Y, b_y, precond, maxit_y, tol_y);
				printf(" conv Y = %s ..",(result)?"false":"true");
			}

			#pragma omp section
			{
				result = Solver::BiCG(A, Z, b_z, precond, maxit_z, tol_z);
				printf(" conv Z = %s ..",(result)?"false":"true");
			}
		}

		//printf("\nSolve time = %d ms\n", (int)solveTime.elapsed());

		for(int i = 0; i < N; i++)
			points[Surface_mesh::Vertex(i)] = Point(X(i), Y(i), Z(i));

		if(isVolumePreservation)
		{
			// Volume preservation
			double new_volume = m->volume();
			m->scale(pow(init_volume / new_volume, 1.0 / 3.0));
		}

		// Undo addition of virtual hole center points
		//untreatBorders(QSurfaceMesh, *U);

	}

	//printf("Smoothing done. (%d ms)\n", (int)timer.elapsed());*/
}

void Smoother::MeanCurvatureFlowExplicit( QSurfaceMesh * m, int numIteration, double step )
{
	//CreateTimer(timer);
	printf("\n\nPerforming Mean Curvature Flow smoothing (Explicit, iterations = %d)...", numIteration);

	for(int iteration = 0; iteration < numIteration; iteration++)
	{
		Surface_mesh::Vertex_property<Point> points = m->vertex_property<Point>("v:point");
		Surface_mesh::Vertex_property<Point> d = m->vertex_property<Point>("v:point_delta");
		Surface_mesh::Vertex_property<Point> disp = m->vertex_property<Point>("v:point_displace");

		Surface_mesh::Vertex_iterator v, vend = m->vertices_end();

		// Original positions, for boundary
		for(v = m->vertices_begin(); v != vend; ++v)
		{
			d[v] = Point(0,0,0);
			disp[v] = Point(0,0,0);
		}

		double init_volume = m->volume();

		double alpha, beta;

		// For each point
		for(v = m->vertices_begin(); v != vend; ++v)
		{
			Point p0 = points[v];

			double totw = 0;

			Surface_mesh::Halfedge_around_vertex_circulator h = m->halfedges(v), hend = m->halfedges(v);
			do{
				Point p1 = points[m->to_vertex(h)];

				double w = 0;

				if(m->is_boundary(h))
					w = 0;
				else
				{
					Vec3d evec = p1 - p0;

					Point p2 = points[m->to_vertex(m->next_halfedge(h))];
					Point p3 = points[m->to_vertex(m->next_halfedge(m->opposite_halfedge(h)))];
					Vec3d d0 = (p0 - p2).normalize();
					Vec3d d1 = (p1 - p2).normalize();
					Vec3d d2 = (p0 - p3).normalize();
					Vec3d d3 = (p1 - p3).normalize();

					double cos_alpha = dot(d0,d1);
					double cos_beta = dot(d2,d3);

					alpha = acos( cos_alpha );
					beta = acos( cos_beta );

					double cot_alpha = cos_alpha / sin(alpha);
					double cot_beta = cos_beta / sin(beta);

					double w = cot_alpha + cot_beta;
					totw += w;
					d[v] += w * evec;
				}

				++h;
			} while( h != hend );

			disp[v] = d[v] * (1.0 / totw);
		}

		// Set vertices to final position
		for(v = m->vertices_begin(); v != vend; ++v)
			points[v] += disp[v];

		// Volume preservation
		double new_volume = m->volume();
		m->scale(pow(init_volume / new_volume, 1.0 / 3.0));
	}

	//printf(" done. (%d ms)\n", (int)timer.elapsed());
}

#if 0
void Smoother::treatBorders(QSurfaceMesh * QSurfaceMesh, Vector<Umbrella> & U)
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

					// Un-flag as border
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

void Smoother::untreatBorders(QSurfaceMesh * QSurfaceMesh, Vector<Umbrella> & U)
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
				mesh->tempUmbrellas[j] = Umbrella(mesh->vd(j));
		}

		U.resize(N);
	}
}

#endif
