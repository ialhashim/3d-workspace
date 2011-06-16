#ifndef PLANE_H
#define PLANE_H

#include "Face.h"

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
	Vec n;
	double d;

	Vec center; // for drawing

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
                        n = from.n.unit();
			d = from.d;

			center = from.center;
		}

		return *this;
	}

	Plane(Face * fromFace)
	{
		CalcPlane(*fromFace->v[0], *fromFace->v[1], *fromFace->v[2]);

                center = projectionOf(Vec(0.01f,0.01f,0.01f)); // recheck
	}

	Plane(const Vec& fromNormal, const Vec& point = Vec())
	{
		n = fromNormal;
		d = -(n * point);

		center = projectionOf(point);
	}

	Plane(const Vec &pta, const Vec &ptb, const Vec &ptc)
	{
		CalcPlane(pta, ptb, ptc);
		center = (pta + ptb + ptc) / 3.0;
	}

	inline void CalcPlane(const Vec &pta, const Vec &ptb, const Vec &ptc)
	{
		n = ((ptb - pta) ^ (ptc - pta)).unit();
		d = -(n * pta);
	}

	Plane inverseDir()
	{
		return Plane(-n, center);
	}

	inline bool IsFront(const Vec &pt) const 
	{
		return (GetPointDistance(pt) > 0.0) ? 1 : 0;
	}

	inline bool IsBack(const Vec& pt) const          
	{
		return !(IsFront(pt));
	}

	inline bool IsOn(const Vec& pt) const
	{
		double dist = GetPointDistance(pt);
		return (dist > -Epsilon && dist < Epsilon) ? 1 : 0;
	}

	double GetPointDistance(const Vec& pt ) const;

	inline bool IsFront(const Vec &v1, const Vec &v2, const Vec &v3) const
	{
		return IsFront(v1) && IsFront(v2) && IsFront(v3);
	}

	inline bool IsBack(const Vec &v1, const Vec &v2, const Vec &v3) const
	{
		return !(IsFront(v1,v2,v3));
	}

	inline bool IsFront(Face * face) const
	{
		return IsFront(face->vec(0), face->vec(1), face->vec(2));
	}

	bool IsInTri(const Vec& p, const Vec &a, const Vec &b, const Vec &c) const;

	Vec projectionOf(const Vec &point);

	void projectLine(Line & line);

        int LineIntersect(const Line& l, Vec & result);
	int LineIntersect(const Vec& start, const Vec& end, const Vec& pointOnPlane, Vec & result);

        int ContourFacet(BaseTriangle * f, Vec & p1, Vec & p2);

	void draw(double extent = 0.25);

	void getTangents(Vec& X, Vec& Y) const;

};

#endif // PLANE_H
