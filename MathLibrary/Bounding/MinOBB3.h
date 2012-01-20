// From Geometric Tools, LLC

#pragma once

#include "SurfaceMesh/Vector.h"
#include <vector>
#include "ConvexHull3.h"
#include "MinOBB2.h"
#include "limits.h"

typedef double					Real;
typedef Vector<Real, 3>			Vector3;
typedef Vector<Real, 2>			Vector2;
#define MAX_REAL std::numeric_limits<Real>::infinity()

// Compute a minimum volume oriented box containing the specified points.

class  MinOBB3
{

public:
	class Box3
	{
	public:
		// Constructor
		Box3();
		Box3( Box3 &box );

		// Regularize
		void normalizeAxis();
		void makeRightHanded();

		// Operator
		bool operator == (Box3& box);
	
		// Proximity
		Vec3d ClosestPoint(const Vec3d& p);
		void ClosestSegment( Box3 other, Vec3d & p, Vec3d & q);

	public:
		Vector3 Center;
		std::vector<Vector3> Axis;
		Vector3 Extent;
	};


public:
    MinOBB3(std::vector<Vector3> &points);
	MinOBB3(Surface_mesh * mesh);

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
	bool isReady;
};

