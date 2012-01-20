#pragma once

#include "Vector.h"

struct HitResult{
	bool hit;
	double distance;

	double u;
	double v;
	int index;

	HitResult(bool isHit = false, double hitDistance = DBL_MAX) : hit(isHit), distance(hitDistance)
	{
		u = -1.0;
		v = -1.0;
		index = -1;
	}

	HitResult(const HitResult& other)
	{
		this->hit = other.hit;
		this->distance = other.distance;
		this->u = other.u;
		this->v = other.v;
		this->index = other.index;
	}

	HitResult& operator= (const HitResult& other)
	{
		this->hit = other.hit;
		this->distance = other.distance;
		this->u = other.u;
		this->v = other.v;
		this->index = other.index;

		return *this;
	}
};

struct Ray
{
	Vec3d origin;
	Vec3d direction;
	int index;

	Ray(const Vec3d & Origin = Vec3d(), const Vec3d & Direction = Vec3d(), int Index = -1) : origin(Origin), index(Index)
	{
		direction = Direction / Direction.norm();
	}

	inline Ray inverse() const { return Ray(origin, -direction); } 

	Ray& operator= (const Ray& other)
	{
		this->origin = other.origin;
		this->direction = other.direction;
		this->index = other.index;

		return *this;
	}
};

/*
* See http://local.wasp.uwa.edu.au/~pbourke/papers/conrec/ for the original
* paper by Paul D. Bourke.
-------------------------------------------------------------------------
Create a contour slice through a 3 vertex facet
Given the normal of the cutting plane "n" and a point on the plane "p0"
Return
0 if the contour plane doesn't cut the facet
2 if it does cut the facet, the contour line segment is p1->p2
-1 for an unexpected occurrence
If a vertex touches the contour plane nothing need to be drawn!?
Note: the following has been written as a "stand along" piece of
code that will work but is far from efficient....
*/
inline int ContourFacet( const Vec3d& planeNormal, const Vec3d& pointOnPlane, 
	const std::vector<Vec3d>& face, Vec3d & p1, Vec3d & p2 )
{
	double A,B,C,D;
	double side[3];

	Vec3d p [3] = {face[0], face[1], face[2]}; 

	/*
	Determine the equation of the plane as
	Ax + By + Cz + D = 0
	*/
	A = planeNormal.x();
	B = planeNormal.y();
	C = planeNormal.z();
	D = -(planeNormal.x()*pointOnPlane.x() + 
		planeNormal.y()*pointOnPlane.y() + 
		planeNormal.z()*pointOnPlane.z());

	/*
	Evaluate the equation of the plane for each vertex
	If side < 0 then it is on the side to be retained
	else it is to be clipped
	*/
	side[0] = A*p[0].x() + B*p[0].y() + C*p[0].z() + D;
	side[1] = A*p[1].x() + B*p[1].y() + C*p[1].z() + D;
	side[2] = A*p[2].x() + B*p[2].y() + C*p[2].z() + D;

	/* Are all the vertices on the same side */
	if (side[0] >= 0 && side[1] >= 0 && side[2] >= 0)	return 0 ;
	if (side[0] <= 0 && side[1] <= 0 && side[2] <= 0)	return 0 ;

	/* Is p0 the only point on a side by itself */
	if ((SIGN(side[0]) != SIGN(side[1])) && (SIGN(side[0]) != SIGN(side[2]))) {
		p1.x() = p[0].x() - side[0] * (p[2].x() - p[0].x()) / (side[2] - side[0]);
		p1.y() = p[0].y() - side[0] * (p[2].y() - p[0].y()) / (side[2] - side[0]);
		p1.z() = p[0].z() - side[0] * (p[2].z() - p[0].z()) / (side[2] - side[0]);
		p2.x() = p[0].x() - side[0] * (p[1].x() - p[0].x()) / (side[1] - side[0]);
		p2.y() = p[0].y() - side[0] * (p[1].y() - p[0].y()) / (side[1] - side[0]);
		p2.z() = p[0].z() - side[0] * (p[1].z() - p[0].z()) / (side[1] - side[0]);
		return(2);
	}

	/* Is p1 the only point on a side by itself */
	if ((SIGN(side[1]) != SIGN(side[0])) && (SIGN(side[1]) != SIGN(side[2]))) {
		p1.x() = p[1].x() - side[1] * (p[2].x() - p[1].x()) / (side[2] - side[1]);
		p1.y() = p[1].y() - side[1] * (p[2].y() - p[1].y()) / (side[2] - side[1]);
		p1.z() = p[1].z() - side[1] * (p[2].z() - p[1].z()) / (side[2] - side[1]);
		p2.x() = p[1].x() - side[1] * (p[0].x() - p[1].x()) / (side[0] - side[1]);
		p2.y() = p[1].y() - side[1] * (p[0].y() - p[1].y()) / (side[0] - side[1]);
		p2.z() = p[1].z() - side[1] * (p[0].z() - p[1].z()) / (side[0] - side[1]);
		return(2);
	}

	/* Is p2 the only point on a side by itself */
	if ((SIGN(side[2]) != SIGN(side[0])) && (SIGN(side[2]) != SIGN(side[1]))) {
		p1.x() = p[2].x() - side[2] * (p[0].x() - p[2].x()) / (side[0] - side[2]);
		p1.y() = p[2].y() - side[2] * (p[0].y() - p[2].y()) / (side[0] - side[2]);
		p1.z() = p[2].z() - side[2] * (p[0].z() - p[2].z()) / (side[0] - side[2]);
		p2.x() = p[2].x() - side[2] * (p[1].x() - p[2].x()) / (side[1] - side[2]);
		p2.y() = p[2].y() - side[2] * (p[1].y() - p[2].y()) / (side[1] - side[2]);
		p2.z() = p[2].z() - side[2] * (p[1].z() - p[2].z()) / (side[1] - side[2]);
		return(2);
	}

	/* Shouldn't get here */
	return 0 ;
}

// Point (P) ------ Triangle (a,b,c)
inline Point ClosestPtPointTriangle(Point p, Point a, Point b, Point c)
{
	// Check if P in vertex region outside A
	Vec3d ab = b - a;
	Vec3d ac = c - a;
	Vec3d ap = p - a;
	double d1 = dot(ab, ap);
	double d2 = dot(ac, ap);
	if (d1 <= 0 && d2 <= 0) return a; // barycentric coordinates (1,0,0)
	// Check if P in vertex region outside B
	Vec3d bp = p - b;
	double d3 = dot(ab, bp);
	double d4 = dot(ac, bp);
	if (d3 >= 0 && d4 <= d3) return b; // barycentric coordinates (0,1,0)
	// Check if P in edge region of AB, if so return projection of P onto AB
	double vc = d1*d4 - d3*d2;
	if (vc <= 0 && d1 >= 0 && d3 <= 0) {
		double v = d1 / (d1 - d3);
		return a + v * ab; // barycentric coordinates (1-v,v,0)
	}
	// Check if P in vertex region outside C
	Vec3d cp = p - c;
	double d5 = dot(ab, cp);
	double d6 = dot(ac, cp);
	if (d6 >= 0 && d5 <= d6) return c; // barycentric coordinates (0,0,1)
	// Check if P in edge region of AC, if so return projection of P onto AC
	double vb = d5*d2 - d1*d6;
	if (vb <= 0 && d2 >= 0 && d6 <= 0) {
		double w = d2 / (d2 - d6);
		return a + w * ac; // barycentric coordinates (1-w,0,w)
	}
	// Check if P in edge region of BC, if so return projection of P onto BC
	double va = d3*d6 - d5*d4;
	if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
		double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		return b + w * (c - b); // barycentric coordinates (0,1-w,w)
	}
	// P inside face region. Compute Q through its barycentric coordinates (u,v,w)
	double denom = 1 / (va + vb + vc);
	double v = vb * denom;
	double w = vc * denom;
	return a + ab * v + ac * w; // = u*a + v*b + w*c, u = va * denom = 1 - v - w
}

inline double distancePointTri(Point p, Point a, Point b, Point c)
{
	return (p - ClosestPtPointTriangle(p,a,b,c)).norm();
}

inline void ClosestPtPointSegment(Point c, Point a, Point b, double &t, Point &d)
{
	Vec3d ab = b - a;
	// Project c onto ab, but deferring divide by Dot(ab, ab)
	t = dot(c - a, ab);
	if (t <= 0.0f) {
		// c projects outside the [a,b] interval, on the a side; clamp to a
		t = 0.0f;
		d = a;
	} else {
		double denom = dot(ab, ab); // Always nonnegative since denom = ||ab||
		if (t >= denom) {
			// c projects outside the [a,b] interval, on the b side; clamp to b
			t = 1.0;
			d = b;
		} else {
			// c projects inside the [a,b] interval; must do deferred divide now
			t = t / denom;
			d = a + t * ab;
		}
	}
}

inline double SqDistPointSegment(Point a, Point b, Point c)
{
	Vec3d ab = b - a, ac = c - a, bc = c - b;
	double e = dot(ac, ab);
	// Handle cases where c projects outside ab
	if (e <= 0.0f) return dot(ac, ac);
	double f = dot(ab, ab);
	if (e >= f) return dot(bc, bc);
	// Handle cases where c projects onto ab
	return dot(ac, ac) - e * e / f;
}

// Clamp n to lie within the range [min, max]
inline double Clamp(double n, double min, double max) {
	if (n < min) return min;
	if (n > max) return max;
	return n;
}

// Computes closest points C1 and C2 of S1(s)=P1+s*(Q1-P1) and
// S2(t)=P2+t*(Q2-P2), returning s and t. Function result is squared
// distance between between S1(s) and S2(t)
inline double ClosestPtSegmentSegment(Point p1, Point q1, Point p2, Point q2,
	double &s, double &t, Point &c1, Point &c2, double EPSILON = 1e-10)
{
	Vec3d d1 = q1 - p1; // Direction Vec3d of segment S1
	Vec3d d2 = q2 - p2; // Direction Vec3d of segment S2
	Vec3d r = p1 - p2;
	double a = dot(d1, d1); // Squared length of segment S1, always nonnegative
	double e = dot(d2, d2); // Squared length of segment S2, always nonnegative
	double f = dot(d2, r);
	// Check if either or both segments degenerate into points
	if (a <= EPSILON && e <= EPSILON) {
		// Both segments degenerate into points
		s = t = 0.0;
		c1 = p1;
		c2 = p2;
		return dot(c1 - c2, c1 - c2);
	}
	if (a <= EPSILON) {
		// First segment degenerates into a point
		s = 0.0;
		t = f / e; // s = 0 => t = (b*s + f) / e = f / e
		t = Clamp(t, 0.0, 1.0);
	} else {
		double c = dot(d1, r);
		if (e <= EPSILON) {
			// Second segment degenerates into a point
			t = 0.0;
			s = Clamp(-c / a, 0.0, 1.0); // t = 0 => s = (b*t - c) / a = -c / a
		} else {
			// The general non-degenerate case starts here
			double b = dot(d1, d2);
			double denom = a*e-b*b; // Always nonnegative
			// If segments not parallel, compute closest point on L1 to L2 and
			// clamp to segment S1. Else pick arbitrary s (here 0)
			if (denom != 0.0) {
				s = Clamp((b*f - c*e) / denom, 0.0, 1.0);
			} else s = 0.0;
			// Compute point on L2 closest to S1(s) using
			// t = dot((P1 + D1*s) - P2,D2) / dot(D2,D2) = (b*s + f) / e
			t = (b*s + f) / e;
			// If t in [0,1] done. Else clamp t, recompute s for the new value
			// of t using s = dot((P2 + D2*t) - P1,D1) / dot(D1,D1)= (t*b - c) / a
			// and clamp s to [0, 1]
			if (t < 0.0) {
				t = 0.0;
				s = Clamp(-c / a, 0.0, 1.0);
			} else if (t > 1.0) {
				t = 1.0;
				s = Clamp((b - c) / a, 0.0, 1.0);
			}
		}
	}
	c1 = p1 + d1 * s;
	c2 = p2 + d2 * t;
	return dot(c1 - c2, c1 - c2);
}
