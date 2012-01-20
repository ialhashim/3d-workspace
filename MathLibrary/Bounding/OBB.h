#pragma once

#include "Surface_mesh.h"

// Eigen: for matrix stuff, eigenvectors
#include "Eigen/Geometry"
#include "Eigen/Eigenvalues"

// For visualization
#include "SimpleDraw.h"

// Based on code by James Gregson: http://jamesgregson.blogspot.com/2011/03/latex-test.html

class OBB {
private:
	Vec3d	m_rot[3];	// rotation matrix for the transformation, stored as rows
	Vec3d	m_pos;		// translation component of the transformation
	Vec3d	m_ext;		// bounding box extents

	Vec3d r,u,f;

	bool isReady;

	// method to set the OBB parameters which produce a box oriented according to
	// the covariance matrix C, which just contains the points pnts
	void build_from_covariance_matrix( Eigen::Matrix3d &C, std::vector<Vec3d> &pnts ){
		// extract the eigenvalues and eigenvectors from C
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(C);
		Eigen::Matrix3d eigvec = es.eigenvectors();

		// find the right, up and forward vectors from the eigenvectors
		r = Vec3d ( eigvec(0,0), eigvec(1,0), eigvec(2,0) );
		u = Vec3d ( eigvec(0,1), eigvec(1,1), eigvec(2,1) );
		f = Vec3d ( eigvec(0,2), eigvec(1,2), eigvec(2,2) );
		r.normalize(); u.normalize(), f.normalize();

		// set the rotation matrix using the eigenvectors
		m_rot[0].x() = r.x(); m_rot[0].y()=u.x(); m_rot[0].z()=f.x();
		m_rot[1].x() = r.y(); m_rot[1].y()=u.y(); m_rot[1].z()=f.y();
		m_rot[2].x() = r.z(); m_rot[2].y()=u.z(); m_rot[2].z()=f.z();

		// now build the bounding box extents in the rotated frame
		double max_dbl = std::numeric_limits<double>::max(), min_dbl = std::numeric_limits<double>::min();
		Vec3d minim(max_dbl, max_dbl, max_dbl);
		Vec3d maxim(min_dbl, min_dbl, min_dbl);

		for( int i = 0; i < (int)pnts.size(); i++ ){
			Vec3d p_prime( dot(r,pnts[i]), dot(u, pnts[i]), dot(f, pnts[i]) );
			minim = minim.minimize( p_prime );
			maxim = maxim.maximize( p_prime );
		}

		// set the center of the OBB to be the average of the 
		// minimum and maximum, and the extents be half of the
		// difference between the minimum and maximum
		Vec3d center = (maxim+minim)*0.5;
		m_pos = Vec3d( dot(m_rot[0],center), dot(m_rot[1], center), dot(m_rot[2],center) );
		m_ext = (maxim-minim)*0.5;

		isReady = true;
	}

public:
	OBB(){ isReady = false; }

	// returns the volume of the OBB, which is a measure of
	// how tight the fit is.  Better OBBs will have smaller 
	// volumes
	double volume(){
		return 8*m_ext[0]*m_ext[1]*m_ext[2];
	}

	// constructs the corner of the aligned bounding box
	// in world space
	void get_bounding_box( Vec3d *p ){
		Vec3d r( m_rot[0][0], m_rot[1][0], m_rot[2][0] );
		Vec3d u( m_rot[0][1], m_rot[1][1], m_rot[2][1] );
		Vec3d f( m_rot[0][2], m_rot[1][2], m_rot[2][2] );
		p[0] = m_pos - r*m_ext[0] - u*m_ext[1] - f*m_ext[2];
		p[1] = m_pos + r*m_ext[0] - u*m_ext[1] - f*m_ext[2];
		p[2] = m_pos + r*m_ext[0] - u*m_ext[1] + f*m_ext[2];
		p[3] = m_pos - r*m_ext[0] - u*m_ext[1] + f*m_ext[2];
		p[4] = m_pos - r*m_ext[0] + u*m_ext[1] - f*m_ext[2];
		p[5] = m_pos + r*m_ext[0] + u*m_ext[1] - f*m_ext[2];
		p[6] = m_pos + r*m_ext[0] + u*m_ext[1] + f*m_ext[2];
		p[7] = m_pos - r*m_ext[0] + u*m_ext[1] + f*m_ext[2];
	}

	// build an OBB from a vector of input points.  This
	// method just forms the covariance matrix and hands
	// it to the build_from_covariance_matrix() method
	// which handles fitting the box to the points
	void build_from_points( std::vector<Vec3d> &pnts ){
		Vec3d mu(0.0, 0.0, 0.0);
		Eigen::Matrix3d C;

		// loop over the points to find the mean point
		// location
		for( int i=0; i<(int)pnts.size(); i++ ){
			mu += pnts[i] / double(pnts.size());
		}

		// loop over the points again to build the 
		// covariance matrix.  Note that we only have
		// to build terms for the upper triangular 
		// portion since the matrix is symmetric
		double cxx=0.0, cxy=0.0, cxz=0.0, cyy=0.0, cyz=0.0, czz=0.0;
		for( int i=0; i < (int)pnts.size(); i++ ){
			Vec3d &p = pnts[i];
			cxx += p.x()*p.x() - mu.x()*mu.x(); 
			cxy += p.x()*p.y() - mu.x()*mu.y(); 
			cxz += p.x()*p.z() - mu.x()*mu.z();
			cyy += p.y()*p.y() - mu.y()*mu.y();
			cyz += p.y()*p.z() - mu.y()*mu.z();
			czz += p.z()*p.z() - mu.z()*mu.z();
		}

		// now build the covariance matrix
		C(0,0) = cxx; C(0,1) = cxy; C(0,2) = cxz;
		C(1,0) = cxy; C(1,1) = cyy; C(1,2) = cyz;
		C(2,0) = cxz; C(2,1) = cyz; C(2,2) = czz;

		// set the OBB parameters from the covariance matrix
		build_from_covariance_matrix( C, pnts );
	}

	// builds an OBB from triangles specified as an array of
	// points with integer indices into the point array. Forms
	// the covariance matrix for the triangles, then uses the
	// method build_from_covariance_matrix() method to fit 
	// the box.  ALL points will be fit in the box, regardless
	// of whether they are indexed by a triangle or not.
	void build_from_triangles( std::vector<Vec3d> &pnts, std::vector<int> &tris ){
		double Ai, Am=0.0;
		Vec3d mu(0.0f, 0.0f, 0.0f), mui;
		Eigen::Matrix3d C;
		double cxx=0.0, cxy=0.0, cxz=0.0, cyy=0.0, cyz=0.0, czz=0.0;

		// loop over the triangles this time to find the
		// mean location
		for( int i=0; i<(int)tris.size(); i+=3 ){
			Vec3d &p = pnts[tris[i+0]];
			Vec3d &q = pnts[tris[i+1]];
			Vec3d &r = pnts[tris[i+2]];
			mui = (p+q+r)/3.0;
			Ai = cross((q-p),(r-p)).norm()/2.0;

			mu += mui*Ai;
			Am += Ai;

			// these bits set the c terms to Am*E[xx], Am*E[xy], Am*E[xz]....
			cxx += ( 9.0*mui.x()*mui.x() + p.x()*p.x() + q.x()*q.x() + r.x()*r.x() )*(Ai/12.0);
			cxy += ( 9.0*mui.x()*mui.y() + p.x()*p.y() + q.x()*q.y() + r.x()*r.y() )*(Ai/12.0);
			cxz += ( 9.0*mui.x()*mui.z() + p.x()*p.z() + q.x()*q.z() + r.x()*r.z() )*(Ai/12.0);
			cyy += ( 9.0*mui.y()*mui.y() + p.y()*p.y() + q.y()*q.y() + r.y()*r.y() )*(Ai/12.0);
			cyz += ( 9.0*mui.y()*mui.z() + p.y()*p.z() + q.y()*q.z() + r.y()*r.z() )*(Ai/12.0);
		}
		// divide out the Am fraction from the average position and 
		// covariance terms
		mu /= Am;
		cxx /= Am; cxy /= Am; cxz /= Am; cyy /= Am; cyz /= Am; czz /= Am;

		// now subtract off the E[x]*E[x], E[x]*E[y], ... terms
		cxx -= mu.x()*mu.x(); cxy -= mu.x()*mu.y(); cxz -= mu.x()*mu.z();
		cyy -= mu.y()*mu.y(); cyz -= mu.y()*mu.z(); czz -= mu.z()*mu.z();

		// now build the covariance matrix
		C(0,0)=cxx; C(0,1)=cxy; C(0,2)=cxz;
		C(1,0)=cxy; C(1,1)=cyy; C(1,2)=cyz;
		C(2,0)=cxz; C(1,2)=cyz; C(2,2)=czz;

		// set the OBB parameters from the covariance matrix
		build_from_covariance_matrix( C, pnts );
	}

	void build_from_mesh( Surface_mesh * mesh )
	{
		// Get points
		std::vector<Vec3d> pnts;
		
		Surface_mesh::Vertex_property<Point> points = mesh->vertex_property<Point>("v:point");
		Surface_mesh::Vertex_iterator vit, vend = mesh->vertices_end();

		for (vit = mesh->vertices_begin(); vit != vend; ++vit)
			pnts.push_back(points[vit]);

		// Get face indices
		std::vector<int> tris;
		Surface_mesh::Face_iterator fit, fend = mesh->faces_end();
		Surface_mesh::Vertex_around_face_circulator fvit;
		Surface_mesh::Vertex v0, v1, v2;

		for(fit = mesh->faces_begin(); fit != fend; ++fit){
			fvit = mesh->vertices(fit);
			v0 = fvit; v1 = ++fvit; v2 = ++fvit;
			tris.push_back(v0.idx()); 
			tris.push_back(v1.idx()); 
			tris.push_back(v2.idx());
		}
		
		build_from_triangles(pnts, tris);
	}

	void draw()
	{
		if(!isReady) return;

		Vec3d pnts[8];
		get_bounding_box(pnts);

		SimpleDraw::IdentifyLine(pnts[0], pnts[1]);
		SimpleDraw::IdentifyLine(pnts[1], pnts[2]);
		SimpleDraw::IdentifyLine(pnts[2], pnts[3]);
		SimpleDraw::IdentifyLine(pnts[3], pnts[0]);

		SimpleDraw::IdentifyLine(pnts[0], pnts[4]);
		SimpleDraw::IdentifyLine(pnts[1], pnts[5]);
		SimpleDraw::IdentifyLine(pnts[2], pnts[6]);
		SimpleDraw::IdentifyLine(pnts[3], pnts[7]);

		SimpleDraw::IdentifyLine(pnts[4], pnts[5]);
		SimpleDraw::IdentifyLine(pnts[5], pnts[6]);
		SimpleDraw::IdentifyLine(pnts[6], pnts[7]);
		SimpleDraw::IdentifyLine(pnts[7], pnts[4]);
	}

	Point center()
	{
		return m_pos;
	}

	Vec3d extents()
	{
		return m_ext;
	}

	std::vector<Vec3d> axis()
	{
		std::vector<Vec3d> result;

		result.push_back(r);
		result.push_back(u);
		result.push_back(f);

		return result;
	}
};
