#include "GCDeformation.h"

// Only needed for one task (when points outside cage case)
#include <Eigen/Geometry>
using namespace Eigen;

GCDeformation::GCDeformation( QSurfaceMesh * forShape, QSurfaceMesh * usingCage )
{
	this->shape = forShape;
	this->cage = usingCage;
	
	orginalCagePos = cage->clonePoints();
	orginalCageNormal = cage->cloneFaceNormals();

	coords.resize(shape->n_vertices());

	// For all points in shape, compute coordinates
	std::vector<Point> shapePoints = shape->clonePoints();
	
	#pragma omp parallel for
	for (int i = 0; i < (int)shapePoints.size(); i++)
	{
		int numRetry = 0;

		GreenCoordiante gc = computeCoordinates( shapePoints[i]);
		
		// Numerical issue, solved by adding small noise to point
		while(!gc.valid){
			double t = 1e-6;
			Point q = shapePoints[i] + Vec3d(uniform(0,t), uniform(0,t), uniform(0,t));
			gc = computeCoordinates(q);

			numRetry++;

			if(numRetry > 10)
				break;
		}

		coords[i].coord_v = gc.coord_v;
		coords[i].coord_n = gc.coord_n;
	}
}

GCDeformation::GreenCoordiante GCDeformation::computeCoordinates(Vec3d point)
{
	GreenCoordiante gc;

	Surface_mesh::Face_iterator fit, fend = cage->faces_end();

	gc.coord_v.resize(cage->n_vertices(), 0);
	gc.coord_n.resize(cage->n_faces(), 0);

	gc.insideCage = true;
	gc.valid = true;

	// For each face in cage
	for(fit = cage->faces_begin(); fit != fend; ++fit)
	{
		Vec3d v[3], s, I, II, N[3], Zero(0,0,0);

		std::vector<Vec3d> facePnts = cage->facePoints(fit);
		std::vector<uint> faceVrts = cage->faceVerts(fit);
		uint fi = Surface_mesh::Face(fit).idx();
		Vec3d n = orginalCageNormal[ fi ];
		
		// 1) First "foreach"
		for (int l = 0; l < 3; l++)
			v[l] = (facePnts[l] - point);

		// 2 ) Assign "p"
		Vec3d p = dot(v[0], n) * n;

		// 3) For each vertex 1, 2, 3
		for (int l = 0; l < 3; l++) 
		{
			int l1 = (l + 1) % 3;

			double DOT = dot((cross((v[l] - p), (v[l1] - p))), n);
			s [l] = DOT < 0.0 ? -1.0 : 1.0;  // Sign

			I [l] = GCTriInt(p,    v[l], v[l1], Zero);
			II[l] = GCTriInt(Zero, v[l1], v[l], Zero);

			N [l] = cross(v[l1] , v[l]).normalized();
		}

		// 4) Psi
		double& psi = gc.coord_n[ fi ];
		for (int k = 0; k < 3; k++)
			psi += s[k] * I[k];

		psi = abs(psi);

		// Robustness check (bug?)
		if(psi != psi)
		{
			gc.valid = false;
			return gc;
		}

		// 5) "w"
		Vec3d w = -psi * n;
		for (int k = 0; k < 3; k++)
			w += II[k] * N[k];

		// 6) Phi
		for (int l = 0; l < 3; l++)
		{
			int l1 = (l + 1) % 3;

			double& phi = gc.coord_v[ faceVrts[l] ];
			phi += dot(N[l1], w) / dot(N[l1], v[l]);
		}
	}

	double coord_v_sum = 0;

	for (int i = 0; i < cage->n_vertices(); ++i)
		coord_v_sum += gc.coord_v[i];

	if(coord_v_sum < 0.5f)	gc.insideCage = false;

	// Check if vertex is exterior to the cage
	if (!gc.insideCage) 
	{
		Surface_mesh::Face closestFace;
		double dist_min = DBL_MAX;

		// find the nearest face (naive, based on Euclid distance)
		for(fit = cage->faces_begin(); fit != fend; ++fit)
		{
			double dist = 0;

			dist += (cage->faceCenter(fit) - point).norm();

			if (dist < dist_min){
				dist_min = dist;
				closestFace = fit;
			}
		}

		// Get face points and vertex indices
		std::vector<Vec3d> facePnts = cage->facePoints(closestFace);
		std::vector<uint> faceVrts = cage->faceVerts(closestFace);
		uint fi = Surface_mesh::Face(closestFace).idx();

		// compute alpha[3] and beta
		Matrix4d A;

		for (int i = 0; i < 3; i++){
			Vector4d fp; fp << facePnts[i].x(), facePnts[i].y(), facePnts[i].z(), 1.0;
			A.col(i) = fp;
		}

		Vector4d fn; fn << orginalCageNormal[ fi ].x(), orginalCageNormal[ fi ].y(), orginalCageNormal[ fi ].z(), 0.0;
		A.col(3) = fn;

		Vector4d pnt; pnt << point.x(), point.y(), point.z(), 1.0;

		Vector4d x = A.fullPivLu().solve(pnt);

		// Set special coordinates
		for (int i = 0; i < 3; i++)
			gc.coord_v[ faceVrts[i] ] += x[i];
		gc.coord_n[ fi ] += x[3];
	}

	return gc;
}

double GCDeformation::GCTriInt(const Vec3d& p, const Vec3d& v1, const Vec3d& v2, const Vec3d& e)
{
	const double alpha = acos(dot((v2 - v1).normalized() , (p - v1).normalized()));
	const double beta  = acos(dot((v1 - p ).normalized() , (v2 -p ).normalized()));
	const double lambda = (p - v1).sqrnorm() * sin(alpha) * sin(alpha);
	const double c      = (p - e).sqrnorm();
	const double theta[2] = { M_PI - alpha, M_PI - alpha - beta };

	Vec2d I(0,0);

	for (int i = 0; i < 2; ++i)
	{
		double S = sin(theta[i]);
		double C = cos(theta[i]);

		double sign = S < 0 ? -1.0 : 0 < S ? 1.0 : 0.0;

		if (sign == 0.0)
			I[i] = 0.0;
		else
		{
			double M = (-sign / 2.0);
			double N = 2 * sqrt(c) * atan((sqrt(c) * C) / sqrt(lambda + (S * S * c)));
			double O = sqrt(lambda);
			double P = (2 * sqrt(lambda) * S * S) / pow(1.0 - C, 2);
			double denom = ( (c*(1+C) + lambda + sqrt((lambda * lambda) + (lambda * c * S * S)) ));
			double Q = (2 * c * C) / denom;
			double R = 1.0 - Q;

			I[i] = M * (N + (O * log(P * R)));
		}
	}

    double ret = (-0.25 / M_PI) * abs(I[0] - I[1] - sqrt(c) * beta);
	
	return ret;
}

void GCDeformation::initDeform()
{
	// Save deformed cage position and face normals
	cage->update_face_normals();
	deformedCageNormal = cage->cloneFaceNormals();
	deformedCagePos = cage->clonePoints();

	// Compute scale factor per face
	Surface_mesh::Face_iterator fit, fend = cage->faces_end();
	for(fit = cage->faces_begin(); fit != fend; ++fit)
	{
		std::vector<uint> faceVrts = cage->faceVerts(fit);
		uint vi0 = faceVrts[0], vi1 = faceVrts[1], vi2 = faceVrts[2];

		const Vec3d& p00 = orginalCagePos[vi0];	const Vec3d& p10 = deformedCagePos[vi0];
		const Vec3d& p01 = orginalCagePos[vi1];	const Vec3d& p11 = deformedCagePos[vi1];
		const Vec3d& p02 = orginalCagePos[vi2];	const Vec3d& p12 = deformedCagePos[vi2];

		Vec3d u0 = p01 - p00; Vec3d u1 = p11 - p10;
		Vec3d v0 = p02 - p00; Vec3d v1 = p12 - p10;

		S.push_back(sqrt(u1.sqrnorm() * v0.sqrnorm() - 2.0f * dot(u1, v1) * dot(u0, v0) + v1.sqrnorm() * u0.sqrnorm()) 
			/ (sqrt(8.0f) * cage->faceArea(fit)));
	}
}

Point GCDeformation::deformedPoint(GreenCoordiante gc)
{
	// Apply deformation
	Vec3d newPoint (0,0,0);

	for (uint i = 0; i < cage->n_vertices(); i++) 
		newPoint += gc.coord_v[i] * deformedCagePos[i];

	for (uint j = 0; j < cage->n_faces(); j++) 
		newPoint += gc.coord_n[j] * S[j] * deformedCageNormal[j];

	return newPoint;
}

void GCDeformation::deform()
{
	initDeform();

	std::vector<Point> shapePoints = shape->clonePoints();
	
	for (uint i = 0; i < shapePoints.size(); i++)
		shapePoints[i] = deformedPoint( coords[i] );

	shape->setFromPoints(shapePoints);
}
