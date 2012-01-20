#pragma once

#include "Surface_mesh.h"
#include "GraphicsLibrary/SpacePartition/Intersection.h"

#define MODIFIED_TRI 5

class BaseTriangle
{
public:
	int index;
	int flag;

	virtual Vec3d vec(int i) const = 0;
	virtual Vec3d normal() const = 0;

	virtual void intersectionTest(const Ray & ray, HitResult & res, bool allowBack = false) const = 0;
};

class Triangle : public BaseTriangle
{
private:
	Vec3d p[3];

public:

	Triangle()
	{
		p[0] = Vec3d();
		p[1] = Vec3d();
		p[2] = Vec3d();

		index = flag = -1;
	}

	Triangle(const Vec3d& point1, const Vec3d& point2, const Vec3d& point3, int tri_index = -1, int tri_flag = -1)
	{
		p[0] = point1;
		p[1] = point2;
		p[2] = point3;

		index = tri_index;
		flag = tri_flag;
	}

	Vec3d vec(int i) const{ return p[i]; }
	Vec3d center() const { return (p[0]+p[1]+p[2]) / 3.0; };

	Vec3d normal() const
	{
		Vec3d n = cross((p[1] - p[0]), (p[2] - p[0]));

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
