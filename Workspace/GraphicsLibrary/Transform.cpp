#include <iostream>

#include "Transform.h"

Transformation Transform3D::findFrom(Vector<Vec> A, Vector<Vec> B)
{
	// Add anchor points, help avoid symmetric issues 
	Vec centerA = Point3D::MidVec(A); Vec a0 = A[0] - centerA, a1 = A[1] - centerA;
	Vec centerB = Point3D::MidVec(B); Vec b0 = B[0] - centerB, b1 = B[1] - centerB;

	Vec anchorA = ((a0 ^ a1).unit() * 2.0f * a0.norm()) + centerA;
	Vec anchorB = ((b0 ^ b1).unit() * 2.0f * b0.norm()) + centerB;

	A.push_back(anchorA);
	B.push_back(anchorB);

	int n = A.size();
	MatrixXd matA(3, n), matB(3, n);

	for(int i = 0; i < n; i++)
	{
		matA(0, i) = A[i].x;	matA(1, i) = A[i].y;	matA(2, i) = A[i].z; // matA (x,y,z)
		matB(0, i) = B[i].x;	matB(1, i) = B[i].y;	matB(2, i) = B[i].z; // matB (x,y,z)
	}

	Matrix4d tm = umeyama(matA, matB, true);	// figure out transformation (from matrix to matrix)
	Transform3d T(tm);							// convert to Eigen::Transofrm3f object

	Quaterniond q;
	Vector3d scaleVec;

	// Rotation & Scale
	extractScaleRotation(T, scaleVec, q);
	qglviewer::Quaternion r(q.x(), q.y(), q.z(), q.w());
	float s = scaleVec(0);

	// Translation
	Vector3d trans(T.translation());
	Vec t = Vec(trans(0), trans(1), trans(2));

	return Transformation(s, r, t);
}

void Transform3D::transformVertices(Mesh * mesh, const Vector<int>& vertices, const Transformation& T)
{
	for(int i = 0; i < (int)vertices.size(); i++)
	{
		Point3D * v = mesh->v(vertices[i]);

		v->set((T.R.rotate(*v * T.S) + T.t));
	}

	mesh->setDirtyVBO(true);
}

void Transform3D::transformVertex(Vec & v, const Transformation& T)
{
	v = T.R.rotate(v * T.S) + T.t;
}

void Transform3D::extractScaleRotation(const Transform3d& tr, Vector3d &scale, Quaterniond &rot)
{
	const double *data = tr.data();

	scale.x() = sqrt(sqr(data[0]) + sqr(data[1]) + sqr(data[2]));
	scale.y() = sqrt(sqr(data[4]) + sqr(data[5]) + sqr(data[6]));
	scale.z() = sqrt(sqr(data[8]) + sqr(data[9]) + sqr(data[10]));

	Transform3d tb = tr;

	double *b = tb.data();

	b[0] = data[0] / scale.x();  b[1] = data[1] / scale.x();   b[2] = data[2] / scale.x();
	b[4] = data[4] / scale.y();  b[5] = data[5] / scale.y();   b[6] = data[6] / scale.y();
	b[8] = data[8] / scale.z();  b[9] = data[9] / scale.z();   b[10] = data[10] / scale.z();

	double t = 1.00 + b[0] + b[5] + b[10];
	if (t > 0.000001) {
		double q = 2.0 * sqrt(t);
		rot.x() = (b[6] - b[9]) / q;
		rot.y() = (b[8] - b[2]) / q;
		rot.z() = (b[1] - b[4]) / q;
		rot.w() = 0.25 * q;
	} else {
		if (b[0] > b[5] && b[0] > b[10]) {
			t = sqrt(1.0 + b[0] - b[5] - b[10]) * 2.0;
			rot.x() = 0.25 * t;
			rot.y() = (b[4] + b[1]) / t;
			rot.z() = (b[2] + b[8]) / t;
			rot.w() = (b[6] - b[9]) / t;
		} else if (b[5] > b[10]) {
			t = sqrt(1.0 + b[5] - b[0] - b[10]) * 2.0;
			rot.y() = 0.25 * t;
			rot.x() = (b[4] + b[1]) / t;
			rot.w() = (b[8] - b[2]) / t;
			rot.z() = (b[6] + b[9]) / t;
		} else {
			t = sqrt(1.0 + b[10] - b[0] - b[5]) * 2.0;
			rot.z() = 0.25 * t;
			rot.w() = (b[1] - b[4]) / t;
			rot.x() = (b[2] + b[8]) / t;
			rot.y() = (b[6] + b[9]) / t;
		}
	}
}
