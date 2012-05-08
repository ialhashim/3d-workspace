#pragma once

#include "Macros.h"
#include "GraphicsLibrary/Mesh/QSegMesh.h"

#include <QVector>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
using namespace Eigen;


class PCA3
{
public:
	PCA3(QVector<Point> points)
	{
		int n = points.size();
		double x, y, z,	xx, yy, zz, xy, xz, yz,
			sumX(0.0), sumY(0.0), sumZ(0.0),
			sumX2(0.0), sumY2(0.0), sumZ2(0.0),
			sumXY(0.0), sumXZ(0.0), sumYZ(0.0);

		foreach(Point p, points)
		{
			x = p.x();
			y = p.y();
			z = p.z();   
			sumX += x / n;
			sumY += y / n;
			sumZ += z / n;
			sumX2 += x * x / n;
			sumY2 += y * y / n;
			sumZ2 += z * z / n;
			sumXY += x * y / n;
			sumXZ += x * z / n;
			sumYZ += y * z / n;
		}
		xx = sumX2 - sumX * sumX;
		yy = sumY2 - sumY * sumY;
		zz = sumZ2 - sumZ * sumZ;
		xy = sumXY - sumX * sumY;
		xz = sumXZ - sumX * sumZ;
		yz = sumYZ - sumY * sumZ;

		Matrix3d covariance;
		covariance << xx, xy, xz,
					  xy, yy, yz,
					  xz, yz, zz;

		// Solve for eigenvalues and eigenvectors.
		// Eigen values are sorted in ascending order,
		SelfAdjointEigenSolver<Eigen::Matrix3d> es(covariance);
		e_values = es.eigenvalues();
		e_vectors = es.eigenvectors();
	}

	Vec3d eigenvalues()
	{
		return E2V(e_values);
	}

	QVector<Vec3d> eigenvectors()
	{
		QVector<Vec3d> vectors;
		for (int i = 2;i >= 0;i--)
		{
			// Descending order
			vectors.push_back(E2V(e_vectors.col(i)));
		}

		return vectors;
	}
	

private:
	Vector3d e_values;
	Matrix3d e_vectors;
};