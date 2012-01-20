#ifndef PLANE_H
#define PLANE_H

#include "Line.h"

#define NO_INTERSECT 0
#define INTERSECT 1
#define ENDPOINT_INTERSECT 2
#define PARALLEL 3

class Plane
{
private:

public:

	// should be private
	Vec3d n;
	double d;
	Vec3d center; // for drawing

	Plane(){d = 0;}

	Plane(const Plane& from)
	{
		n = from.n;
		d = from.d;

		center = from.center;
	}

	Plane& operator= (const Plane& from)
	{  
		if (this != &from) 
		{
			n = from.n;
			d = from.d;

			center = from.center;
		}

		return *this;
	}
	
	Plane(const Vec3d& fromNormal, const Vec3d& point = Vec3d(0,0,0))
	{
		n = fromNormal.normalized();
		d = -dot(n, point);
		center = point;

		if(point.x() == 0 && point.x() == point.y() && point.y() == point.z())
		{
			d = 0;
			center = Vec3d(0,0,0);
		}
	}

	Plane(const Vec3d &pta, const Vec3d &ptb, const Vec3d &ptc)
	{
		CalcPlane(pta, ptb, ptc);
		center = (pta + ptb + ptc) / 3.0;
	}

	inline void CalcPlane(const Vec3d &pta, const Vec3d &ptb, const Vec3d &ptc)
	{
		n = cross((ptb - pta), (ptc - pta)).normalized();
		d = -dot(n, pta);
	}

	Plane inverseDir()
	{
		return Plane(-n, center);
	}

	inline bool IsFront(const Vec3d &pt) const 
	{
		return (GetPointDistance(pt) > 0.0) ? 1 : 0;
	}

	inline bool IsBack(const Vec3d& pt) const          
	{
		return !(IsFront(pt));
	}

	inline bool IsOn(const Vec3d& pt, double Epsilon = 1e-10) const
	{
		double dist = GetPointDistance(pt);
		return (dist > -Epsilon && dist < Epsilon) ? 1 : 0;
	}

	double GetPointDistance(const Vec3d& pt ) const;

	inline bool IsFront(const Vec3d &v1, const Vec3d &v2, const Vec3d &v3) const
	{
		return IsFront(v1) && IsFront(v2) && IsFront(v3);
	}

	inline bool IsBack(const Vec3d &v1, const Vec3d &v2, const Vec3d &v3) const
	{
		return !(IsFront(v1,v2,v3));
	}
	
	bool IsInTri(const Vec3d& p, const Vec3d &a, const Vec3d &b, const Vec3d &c) const;

	Vec3d projectionOf(const Vec3d &point);

	void projectLine(Line & line);

	int LineIntersect(const Line& l, Vec3d & result = Vec3d());
	int LineIntersect(const Vec3d& start, const Vec3d& end, const Vec3d& pointOnPlane, Vec3d & result, double Epsilon = 1e-10 );
	int ContourFacet(Vec3d a, Vec3d b, Vec3d c, Vec3d & p1, Vec3d & p2);

	void draw(double extent = 0.25);

	void getTangents(Vec3d& X, Vec3d& Y) const;

	static bool isSameSide(const Vec3d& p1, const Vec3d& p2, const Vec3d& a, const Vec3d& b)
	{
		Vec3d cp1 = cross((b-a) , (p1-a));
		Vec3d cp2 = cross((b-a) , (p2-a));

		return dot(cp1, cp2) >= 0;
	}

	static Vec3d orthogonalVector(const Vec3d& n) {
		if ((abs(n.y()) >= 0.9f * abs(n.x())) && 
			abs(n.z()) >= 0.9f * abs(n.x())) return Vec3d(0.0f, -n.z(), n.y());
		else if ( abs(n.x()) >= 0.9f * abs(n.y()) && 
			abs(n.z()) >= 0.9f * abs(n.y()) ) return Vec3d(-n.z(), 0.0f, n.x());
		else return Vec3d(-n.y(), n.x(), 0.0f);
	}

	Vec3d reflection(const Vec3d& v);
};

#endif // PLANE_H
