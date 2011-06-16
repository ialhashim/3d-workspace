#pragma once

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
	Vec origin;
	Vec direction;
	int index;

	Ray(const Vec & Origin = Vec(), const Vec & Direction = Vec(), int Index = -1) : origin(Origin), index(Index)
	{
		direction = Direction.unit();
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
