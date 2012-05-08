#pragma once

#include "Surface_mesh.h"

// Special math library
//#include "OBB2_FloatMath.h"
//#include "OBB2_FloatMath.inl"
#include "OBB_Volume_math.h"

// Eigen: for rotations
#include "Eigen/Geometry"
using namespace Eigen;

class OBB_Volume{

private:
	// computes the OBB for this set of points relative to this transform matrix.
	void computeOBB(size_t vcount,const REAL *points,size_t pstride,REAL *sides,REAL *matrix)
	{
		const char *src = (const char *) points;

		double max_dbl = std::numeric_limits<double>::max(), min_dbl = std::numeric_limits<double>::min();
		REAL bmin[3] = { max_dbl, max_dbl, max_dbl };
		REAL bmax[3] = { min_dbl, min_dbl, min_dbl };

		for (size_t i=0; i<vcount; i++)
		{
			const REAL *p = (const REAL *) src;
			REAL t[3];

			fm_inverseRT(matrix, p, t ); // inverse rotate translate

			if ( t[0] < bmin[0] ) bmin[0] = t[0];
			if ( t[1] < bmin[1] ) bmin[1] = t[1];
			if ( t[2] < bmin[2] ) bmin[2] = t[2];

			if ( t[0] > bmax[0] ) bmax[0] = t[0];
			if ( t[1] > bmax[1] ) bmax[1] = t[1];
			if ( t[2] > bmax[2] ) bmax[2] = t[2];

			src+=pstride;
		}

		REAL center[3];

		sides[0] = bmax[0]-bmin[0];
		sides[1] = bmax[1]-bmin[1];
		sides[2] = bmax[2]-bmin[2];

		center[0] = sides[0]*0.5f+bmin[0];
		center[1] = sides[1]*0.5f+bmin[1];
		center[2] = sides[2]*0.5f+bmin[2];

		REAL ocenter[3];

		fm_rotate(matrix,center,ocenter);

		matrix[12]+=ocenter[0];
		matrix[13]+=ocenter[1];
		matrix[14]+=ocenter[2];

	}

	void fm_computeBestFitOBB(size_t vcount,const REAL *points,size_t pstride,REAL *sides,REAL *matrix,bool bruteForce)
	{
		REAL plane[4];
		fm_computeBestFitPlane(vcount,points,pstride,0,0,plane);
		fm_planeToMatrix(plane,matrix);
		computeOBB( vcount, points, pstride, sides, matrix );

		REAL refmatrix[16];
		memcpy(refmatrix,matrix,16*sizeof(REAL));

		REAL volume = sides[0]*sides[1]*sides[2];
		if ( bruteForce )
		{
			for (REAL a=10; a<180; a+=10)
			{
				REAL quat[4];
				fm_eulerToQuat(0,a*FM_DEG_TO_RAD,0,quat);
				REAL temp[16];
				REAL pmatrix[16];
				fm_quatToMatrix(quat,temp);
				fm_matrixMultiply(temp,refmatrix,pmatrix);
				REAL psides[3];
				computeOBB( vcount, points, pstride, psides, pmatrix );
				REAL v = psides[0]*psides[1]*psides[2];
				if ( v < volume )
				{
					volume = v;
					memcpy(matrix,pmatrix,sizeof(REAL)*16);
					sides[0] = psides[0];
					sides[1] = psides[1];
					sides[2] = psides[2];
				}
			}
		}
	}

	void fm_computeBestFitOBB(size_t vcount,const REAL *points,size_t pstride,REAL *sides,REAL *pos,REAL *quat,bool bruteForce)
	{
		REAL matrix[16];
		fm_computeBestFitOBB(vcount,points,pstride,sides,matrix,bruteForce);
		fm_getTranslation(matrix,pos);
		fm_matrixToQuat(matrix,quat);
	}

public:

	Vec3d sides;
	Vec3d translation;
	Vec4d rotation; // as quat.

	bool isReady;

	OBB_Volume(Surface_mesh * mesh)
	{
		// Get points
		std::vector<Vec3d> pnts;

		Surface_mesh::Vertex_property<Point> points = mesh->vertex_property<Point>("v:point");
		Surface_mesh::Vertex_iterator vit, vend = mesh->vertices_end();

		for (vit = mesh->vertices_begin(); vit != vend; ++vit)
			pnts.push_back(points[vit]);

		sides = translation = Vec3d(0,0,0);
		rotation = Vec4d(0,0,0,0);

		fm_computeBestFitOBB(pnts.size(), &pnts.front()[0], sizeof(Vec3d), &sides[0], &translation[0], &rotation[0], true);

		isReady = true;
	}

	void draw()
	{
		if(!isReady) return;

		Vector3f p(translation.x(), translation.y(), translation.z());
		Vector4f r(rotation[0], rotation[1], rotation[2], rotation[3]);

		Transform<float,3,Affine> t = Translation3f(p) * Quaternionf(r);

		glPushMatrix();
		glMultMatrixf(t.data());

		SimpleDraw::DrawBox(Vec3d(0,0,0), sides.x()/2, sides.y()/2, sides.z()/2);

		glPopMatrix();
	}

	std::vector<Vec3d> axis()
	{
		std::vector<Vec3d> result;

		Vector3f p(0, 0, 0);
		Vector4f r(rotation[0], rotation[1], rotation[2], rotation[3]);

		Transform<float,3,Affine> t = Translation3f(p) * Quaternionf(r);

		Vector3f xAxis(1,0,0);	xAxis = t * xAxis;
		Vector3f yAxis(0,1,0);	yAxis = t * yAxis;
		Vector3f zAxis(0,0,1);	zAxis = t * zAxis;

		result.push_back(Vec3d(xAxis[0], xAxis[1], xAxis[2]));
		result.push_back(Vec3d(yAxis[0], yAxis[1], yAxis[2]));
		result.push_back(Vec3d(zAxis[0], zAxis[1], zAxis[2]));

		return result;
	}

	Point center()
	{
		return translation;
	}

	Vec3d extents()
	{
		return sides * 0.5;
	}
};

