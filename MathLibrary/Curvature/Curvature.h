/* 
	method taken from: 
	Princeton University file 'TriMesh_curvature.cc'
	Rusinkiewicz, Szymon. "Estimating Curvatures and Their Derivatives on Triangle Meshes," 2004
*/

#pragma once

#include "QSurfaceMesh.h"
typedef Surface_mesh::Vertex Vertex;

class Curvature
{

private:
	std::vector<double> curv1, curv2;
	std::vector<Vec4d> dcurv;
	std::vector<Vec3d> pdir1, pdir2;
	std::vector<Vec3d> cornerareas;
	std::vector<double> pointareas;

public:
	// from 'TriMesh2' by "Szymon Rusinkiewicz" - Finite-differences approach
	void computePrincipalCurvatures(Surface_mesh * src_mesh);

	// derivatives of curvature
	void computeDerivativesCurvatures(Surface_mesh * src_mesh);

	// Compute per-vertex point areas
	void computePointAreas(Surface_mesh * src_mesh);

	void rot_coord_sys(const Point &old_u, const Point &old_v, const Point &new_norm, Point &new_u, Point &new_v);

	void proj_curv(const Point &old_u, const Point &old_v, double old_ku, double old_kuv, double old_kv, 
					const Point &new_u, const Point &new_v, double &new_ku, double &new_kuv, double &new_kv);

	void proj_dcurv(const Point &old_u, const Point &old_v, const Vec4d old_dcurv, 
					const Point &new_u, const Point &new_v, Vec4d & new_dcurv);

	void diagonalize_curv(const Point &old_u, const Point &old_v, double ku, double kuv, double kv, 
							const Point &new_norm, Point &pdir1, Point &pdir2, double &k1, double &k2);
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

// Stuff
#define SWAP(x, y, T) do { T temp##x##y = x; x = y; y = temp##x##y; } while (0)
