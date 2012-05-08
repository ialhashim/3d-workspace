// From Geometric Tools, LLC

#pragma once

#include <vector>
#include "ConvexHull3.h"
#include "MinOBB2.h"
#include "Box3.h"

typedef double					Real;
typedef Vector<Real, 3>			Vector3;
typedef Vector<Real, 2>			Vector2;
#define MAX_REAL DOUBLE_INFINITY

// Compute a minimum volume oriented box containing the specified points.

class  MinOBB3
{
public:
    MinOBB3(std::vector<Vector3> &points);
	MinOBB3(QSurfaceMesh * mesh);

	void computeMinOBB( std::vector<Vector3> &points );
	void GenerateComplementBasis (Vector3& u, Vector3& v, const Vector3& w);
	void getCorners( std::vector<Vector3> &pnts );

	Vec3d ClosestPtPointOBB(const Vec3d& p);

private:
	class EdgeKey
	{
	public:
		EdgeKey (int v0 = -1, int v1 = -1);

		bool operator< (const EdgeKey& key) const;
		operator size_t () const;

		int V[2];
	};
    

public:
	Box3 mMinBox;
	Surface_mesh * mesh;
	bool isReady;
};

