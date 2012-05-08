#pragma once

const double TOLERANCE = 0.0000000100;

#include <Eigen/Geometry>
using namespace Eigen;

class DualQuat
{
public:
	DualQuat(void)
	{
		ClearZero();
	}

	// input: unit quaternion 'q', translation vector 't'
	DualQuat(const Quaterniond &qt, const Vector3d &t)
	{
		Vector4d q(qt.w(), qt.x(), qt.y(), qt.z());

		// non-dual part (just copy q):
		for (int i=0; i<4; i++) d[0][i] = q[i];

		// dual part:
		d[1][0] = -0.5*(t[0]*q[1] + t[1]*q[2] + t[2]*q[3]);
		d[1][1] = 0.5*( t[0]*q[0] + t[1]*q[3] - t[2]*q[2]);
		d[1][2] = 0.5*(-t[0]*q[3] + t[1]*q[0] + t[2]*q[1]);
		d[1][3] = 0.5*( t[0]*q[2] - t[1]*q[1] + t[2]*q[0]);
	}
	~DualQuat(void) {}

	void ClearZero(void)
	{
		for (int i=0; i<4; i++) {
			d[0][i] = 0.0;
			d[1][i] = 0.0;
		}
	}

	// Binary operators
	friend DualQuat operator+(const DualQuat& u, const DualQuat& v)
	{
		DualQuat dq;
		for (int i=0; i<4; i++)
		{
			dq.d[0][i] = u.d[0][i] + v.d[0][i];
			dq.d[1][i] = u.d[1][i] + v.d[1][i];
		}
		return dq;
	}
	friend DualQuat operator-(const DualQuat& u, const DualQuat& v)
	{
		DualQuat dq;
		for (int i=0; i<4; i++)
		{
			dq.d[0][i] = u.d[0][i] - v.d[0][i];
			dq.d[1][i] = u.d[1][i] - v.d[1][i];
		}
		return dq;
	}
	friend DualQuat operator*(const double s, const DualQuat& u)
	{
		DualQuat dq;
		for (int i=0; i<4; i++)
		{
			dq.d[0][i] = u.d[0][i] * s;
			dq.d[1][i] = u.d[1][i] * s;
		}
		return dq;
	}
	friend DualQuat operator*(const DualQuat& u, const double s)
	{
		return s * u;
	}
	DualQuat& operator+=(const DualQuat& rd)
	{
		for (int i=0; i<4; i++)
		{
			d[0][i] += rd.d[0][i];
			d[1][i] += rd.d[1][i];
		}
		return *this;
	}
	DualQuat& operator-=(const DualQuat& rd)
	{
		for (int i=0; i<4; i++)
		{
			d[0][i] -= rd.d[0][i];
			d[1][i] -= rd.d[1][i];
		}
		return *this;
	}
	DualQuat& operator/=(const double s)
	{
		for (int i=0; i<4; i++)
		{
			d[0][i] /= s;
			d[1][i] /= s;
		}
		return *this;
	}

	void GetQuat(Quaterniond &qt, Vector3d &t)
	{
		Vector4d q;

		// regular quaternion (just copy the non-dual part):
		for (int i=0; i<4; i++) q[i] = d[0][i];

		qt.w() = q[0];
		qt.x() = q[1];
		qt.y() = q[2];
		qt.z() = q[3];

		// translation vector:
		t[0] = 2.0*(-d[1][0]*d[0][1] + d[1][1]*d[0][0] - d[1][2]*d[0][3] + d[1][3]*d[0][2]);
		t[1] = 2.0*(-d[1][0]*d[0][2] + d[1][1]*d[0][3] + d[1][2]*d[0][0] - d[1][3]*d[0][1]);
		t[2] = 2.0*(-d[1][0]*d[0][3] - d[1][1]*d[0][2] + d[1][2]*d[0][1] + d[1][3]*d[0][0]);
	}

	void SetTransform(const Matrix3d &Rt, const Vector3d &t)
	{
		const double * R = Rt.data();

		// non-dual part (compute quaternion):
		double trace = R[0] + R[4] + R[8] + 1.0;
		if (trace > TOLERANCE) {
			double s = 0.5 / sqrt(trace);
			d[0][0] = 0.25 / s;						// w0
			d[0][1] = (R[7] - R[5]) * s;		// x0
			d[0][2] = (R[2] - R[6]) * s;		// y0
			d[0][3] = (R[3] - R[1]) * s;		// z0
		} else {
			if (R[0]>R[4] && R[0]>R[8]) {
				double s = 0.5 / sqrt(1.0 + R[0] - R[4] - R[8]);
				d[0][0] = (R[5] - R[7]) * s;	// w0
				d[0][1] = 0.25 / s;					// x0
				d[0][2] = (R[1] + R[3]) * s;	// y0
				d[0][3] = (R[2] + R[6]) * s;	// z0
			} else if (R[4] > R[8]) {
				double s = 0.5 / sqrt(1.0 + R[4] - R[0] - R[8]);
				d[0][0] = (R[2] - R[6]) * s;	// w0
				d[0][1] = (R[1] + R[3]) * s;	// x0
				d[0][2] = 0.25 / s;					// y0
				d[0][3] = (R[5] + R[7]) * s;	// z0
			} else {
				double s = 0.5 / sqrt(1.0 + R[8] - R[0] - R[4]);
				d[0][0] = (R[1] - R[3]) * s;	// w0
				d[0][1] = (R[2] + R[6]) * s;	// x0
				d[0][2] = (R[5] + R[7]) * s;	// y0
				d[0][3] = 0.25 / s;					// z0
			}
		}
		// dual part:
		d[1][0] = -0.5*(t[0]*d[0][1] + t[1]*d[0][2] + t[2]*d[0][3]);	// we
		d[1][1] = 0.5*( t[0]*d[0][0] + t[1]*d[0][3] - t[2]*d[0][2]);	// xe
		d[1][2] = 0.5*(-t[0]*d[0][3] + t[1]*d[0][0] + t[2]*d[0][1]);	// ye
		d[1][3] = 0.5*( t[0]*d[0][2] - t[1]*d[0][1] + t[2]*d[0][0]);	// ze
	}

	void SetRotation(const Matrix3d &Rt)
	{
		const double * R = Rt.data();

		// non-dual part (compute quaternion):
		double trace = R[0] + R[4] + R[8] + 1.0;
		if (trace > TOLERANCE) {
			double s = 0.5 / sqrt(trace);
			d[0][0] = 0.25 / s;						// w0
			d[0][1] = (R[7] - R[5]) * s;		// x0
			d[0][2] = (R[2] - R[6]) * s;		// y0
			d[0][3] = (R[3] - R[1]) * s;		// z0
		} else {
			if (R[0]>R[4] && R[0]>R[8]) {
				double s = 0.5 / sqrt(1.0 + R[0] - R[4] - R[8]);
				d[0][0] = (R[5] - R[7]) * s;	// w0
				d[0][1] = 0.25 / s;					// x0
				d[0][2] = (R[1] + R[3]) * s;	// y0
				d[0][3] = (R[2] + R[6]) * s;	// z0
			} else if (R[4] > R[8]) {
				double s = 0.5 / sqrt(1.0 + R[4] - R[0] - R[8]);
				d[0][0] = (R[2] - R[6]) * s;	// w0
				d[0][1] = (R[1] + R[3]) * s;	// x0
				d[0][2] = 0.25 / s;					// y0
				d[0][3] = (R[5] + R[7]) * s;	// z0
			} else {
				double s = 0.5 / sqrt(1.0 + R[8] - R[0] - R[4]);
				d[0][0] = (R[1] - R[3]) * s;	// w0
				d[0][1] = (R[2] + R[6]) * s;	// x0
				d[0][2] = (R[5] + R[7]) * s;	// y0
				d[0][3] = 0.25 / s;					// z0
			}
		}
		// dual part:
		d[1][0] = d[1][1] = d[1][2] = d[1][3] = 0.0;
	}

	void GetTransform(Matrix3d &Rt, Vector3d &t)
	{
		double * R = Rt.data();

		// 3x3 rotation matrix:
		double Nq = (d[0][0]*d[0][0] + d[0][1]*d[0][1] + d[0][2]*d[0][2] + d[0][3]*d[0][3]);
		double c  = (Nq > 0.0) ? (2.0 / Nq) : 0.0;
		double xc = d[0][1]*c,		yc = d[0][2]*c,		zc = d[0][3]*c;
		double wx = d[0][0]*xc,		wy = d[0][0]*yc,	wz = d[0][0]*zc;
		double xx = d[0][1]*xc,		xy = d[0][1]*yc,	xz = d[0][1]*zc;
		double yy = d[0][2]*yc,		yz = d[0][2]*zc,	zz = d[0][3]*zc;
		R[0]	= 1.0 - yy - zz;
		R[1]	= xy - wz;
		R[2]	= xz + wy;
		R[3]	= xy + wz;
		R[4]	= 1.0 - xx - zz;
		R[5]	= yz - wx;
		R[6]	= xz - wy;
		R[7]	= yz + wx;
		R[8]	= 1.0 - xx - yy;
		// translation vector:
		double wc = d[0][0]*c;
		t[0] = -d[1][0]*xc + d[1][1]*wc - d[1][2]*zc + d[1][3]*yc;
		t[1] = -d[1][0]*yc + d[1][1]*zc + d[1][2]*wc - d[1][3]*xc;
		t[2] = -d[1][0]*zc - d[1][1]*yc + d[1][2]*xc + d[1][3]*wc;
	}

	void GetRotation(Matrix3d &Rt)
	{
		double * R = Rt.data();

		// 3x3 rotation matrix:
		double Nq = (d[0][0]*d[0][0] + d[0][1]*d[0][1] + d[0][2]*d[0][2] + d[0][3]*d[0][3]);
		double c  = (Nq > 0.0) ? (2.0 / Nq) : 0.0;
		double xc = d[0][1]*c,		yc = d[0][2]*c,		zc = d[0][3]*c;
		double wx = d[0][0]*xc,		wy = d[0][0]*yc,	wz = d[0][0]*zc;
		double xx = d[0][1]*xc,		xy = d[0][1]*yc,	xz = d[0][1]*zc;
		double yy = d[0][2]*yc,		yz = d[0][2]*zc,	zz = d[0][3]*zc;
		R[0]	= 1.0 - yy - zz;
		R[1]	= xy - wz;
		R[2]	= xz + wy;
		R[3]	= xy + wz;
		R[4]	= 1.0 - xx - zz;
		R[5]	= yz - wx;
		R[6]	= xz - wy;
		R[7]	= yz + wx;
		R[8]	= 1.0 - xx - yy;
	}

	void GetMatrix(Matrix4d &M)
	{
		double * Mat = M.data();

		// 3x3 rotation matrix:
		double Nq = (d[0][0]*d[0][0] + d[0][1]*d[0][1] + d[0][2]*d[0][2] + d[0][3]*d[0][3]);
		double c  = (Nq > 0.0) ? (2.0 / Nq) : 0.0;
		double xc = d[0][1]*c,		yc = d[0][2]*c,		zc = d[0][3]*c;
		double sx = d[0][0]*xc,		sy = d[0][0]*yc,	sz = d[0][0]*zc;
		double xx = d[0][1]*xc,		xy = d[0][1]*yc,	xz = d[0][1]*zc;
		double yy = d[0][2]*yc,		yz = d[0][2]*zc,	zz = d[0][3]*zc;
		Mat[0]	= 1.0 - (yy + zz);
		Mat[1]	= xy - sz;
		Mat[2]	= xz + sy;
		Mat[4]	= xy + sz;
		Mat[5]	= 1.0 - (xx + zz);
		Mat[6]	= yz - sx;
		Mat[8]	= xz - sy;
		Mat[9]	= yz + sx;
		Mat[10]	= 1.0 - (xx + yy);
		// translation vector:
		double wc = d[0][0]*c;
		Mat[3] = -d[1][0]*xc + d[1][1]*wc - d[1][2]*zc + d[1][3]*yc;
		Mat[7] = -d[1][0]*yc + d[1][1]*zc + d[1][2]*wc - d[1][3]*xc;
		Mat[11] = -d[1][0]*zc - d[1][1]*yc + d[1][2]*xc + d[1][3]*wc;
		Mat[12]	= Mat[13] = Mat[14] = 0.0;
		Mat[15]	= 1.0;
	}

private:
	double d[2][4];	// data
};
