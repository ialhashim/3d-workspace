#pragma once

#include "Mesh.h"

struct Transformation
{
	float S;	// scale
        qglviewer::Quaternion R; // rotation
	Vec t;		// translation

        Transformation(float scale, const qglviewer::Quaternion& rotation, const Vec& translation)
	{
		S = scale;
		R = rotation;
		t = translation;
	}
};

class Transform3D
{
public:
	
	typedef Transform<double, 3, Affine> Transform3d;

	// Eigen library realated
	static void extractScaleRotation(const Transform3d& tr, Vector3d &scale, Quaterniond &rot);
	inline static float sqr(float a){return a * a;}

	// PCA and SVD
	static Transformation findFrom(Vector<Vec> A, Vector<Vec> B);
	static void transformVertices(Mesh * mesh, const Vector<int>& vertices, const Transformation& T);
	static void transformVertex(Vec & v, const Transformation& T);
};

