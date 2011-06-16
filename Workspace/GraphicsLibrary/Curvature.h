/* 
	Second method taken from: 
	Princeton University file 'TriMesh_curvature.cc'
	Rusinkiewicz, Szymon. "Estimating Curvatures and Their Derivatives on Triangle Meshes," 2004
*/

#pragma once

#include "Mesh.h"
#include "Vec4.h"

class Curvature
{
private:
	double AREA(double pp1[3], double pp2[3], double pp3[3]);

public:
	// 1st method: from 'MeshViewer' by "Yutaka Ohtake" - Uses one ring neighborhood
	void computePrincipalCurvatures2(Mesh * mesh);

	// 2nd method: from 'TriMesh2' by "Szymon Rusinkiewicz" - Finite-differences approach
	void computePrincipalCurvatures(Mesh * src_mesh);

	// derivatives of curvature
	void computeDerivativesCurvatures(Mesh * src_mesh);

	static void rot_coord_sys(const Vec &old_u, const Vec &old_v, const Vec &new_norm, Vec &new_u, Vec &new_v);

	void proj_curv(const Vec &old_u, const Vec &old_v, double old_ku, double old_kuv, double old_kv, 
					const Vec &new_u, const Vec &new_v, double &new_ku, double &new_kuv, double &new_kv);

	void proj_dcurv(const Vec &old_u, const Vec &old_v, const Vec4 old_dcurv, 
					const Vec &new_u, const Vec &new_v, Vec4 & new_dcurv);

	void diagonalize_curv(const Vec &old_u, const Vec &old_v, double ku, double kuv, double kv, 
							const Vec &new_norm, Vec &pdir1, Vec &pdir2, double &k1, double &k2);

	// Compute per-vertex point areas
	void computePointAreas(Mesh * src_mesh);

public:
	Mesh * mesh;

	// STORED DATA - 1st method
	Vector<double> k_max,k_min;

	// STORED DATA - 2nd method
	Vector<double> curv1, curv2;
	Vector<Vec4> dcurv;
	Vector<Vec> cornerareas;
	Vector<double> pointareas;

	// STORED DATA - shared
	Vector<Vec> pdir1, pdir2;

	// Make copies of mesh (not efficient!)
	Vector<Point3D> vertices;
	Vector<Normal> normals;
	HashMap<int, Face* > faces;

	// Modifiy
	void smoothPrincipalCurvatures(int iterations = 1.0);

};

// Macros needed for Princeton code:

// i+1 and i-1 modulo 3
// This way of computing it tends to be faster than using %
#define NEXT_Index(i) ((i)<2 ? (i)+1 : (i)-2)
#define PREV_Index(i) ((i)>0 ? (i)-1 : (i)+2)

// Let gcc optimize conditional branches a bit better...
#ifndef likely
#  if !defined(__GNUC__) || (__GNUC__ == 2 && __GNUC_MINOR__ < 96)
#    define likely(x) (x)
#    define unlikely(x) (x)
#  else
#    define likely(x)   (__builtin_expect((x), 1))
#    define unlikely(x) (__builtin_expect((x), 0))
#  endif
#endif

// Perform LDL^T decomposition of a symmetric positive definite matrix.
// Like Cholesky, but no square roots.  Overwrites lower triangle of matrix.
template <class T, int N>
static inline bool ldltdc(T A[N][N], T rdiag[N])
{
	T v[N-1];
	for (int i = 0; i < N; i++) {
		for (int k = 0; k < i; k++)
			v[k] = A[i][k] * rdiag[k];
		for (int j = i; j < N; j++) {
			T sum = A[i][j];
			for (int k = 0; k < i; k++)
				sum -= v[k] * A[j][k];
			if (i == j) {
				if (unlikely(sum <= T(0)))
					return false;
				rdiag[i] = T(1) / sum;
			} else {
				A[j][i] = sum;
			}
		}
	}

	return true;
}

// Solve Ax=B after ldltdc
template <class T, int N>
static inline void ldltsl(T A[N][N], T rdiag[N], T B[N], T x[N])
{
	int i;
	for (i = 0; i < N; i++) {
		T sum = B[i];
		for (int k = 0; k < i; k++)
			sum -= A[i][k] * x[k];
		x[i] = sum * rdiag[i];
	}
	for (i = N - 1; i >= 0; i--) {
		T sum = 0;
		for (int k = i + 1; k < N; k++)
			sum += A[k][i] * x[k];
		x[i] -= sum * rdiag[i];
	}
}
