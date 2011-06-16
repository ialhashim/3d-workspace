#pragma once

#include "Point.h"
#include "Intersection.h"

#define MODIFIED_TRI 5

class BaseTriangle
{
public:
	int index;
	int flag;

	virtual Vec vec(int i) const = 0;
	virtual Vec normal() const = 0;

	virtual void intersectionTest(const Ray & ray, HitResult & res, bool allowBack = false) const = 0;
};

class Triangle : public BaseTriangle
{
private:
	Vec p[3];

public:

	Triangle()
	{
		p[0] = Vec();
		p[1] = Vec();
		p[2] = Vec();

		index = flag = -1;
	}

	Triangle(const Vec& point1, const Vec& point2, const Vec& point3, int tri_index = -1, int tri_flag = -1)
	{
		p[0] = point1;
		p[1] = point2;
		p[2] = point3;

		index = tri_index;
		flag = tri_flag;
	}

	Vec vec(int i) const{ return p[i]; }
	Vec center() const { return (p[0]+p[1]+p[2]) / 3.0; };

	Vec normal() const
	{
		Vec n = (p[1] - p[0]) ^ (p[2] - p[0]);

		double length = n.norm();

		if(length < 1.0E-10)
			return n;
		else
			return n /= length;
	}

	inline float edgeLenY(int e) const;

	void unloop(float threshold, Triangle & other);

	Triangle shiftY(float offsetY);

	void intersectionTest(const Ray & ray, HitResult & res, bool allowBack = false) const;
};
