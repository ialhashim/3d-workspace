#ifndef CIRCLE_H
#define CIRCLE_H

#include "Point.h"
#include "Color4.h"

class Circle
{
private:
	float radius;

	int numSides;

	Vec normal;
	Vec center;

	Vector<Vec> point;

public:
	Circle();

	Circle(float from_radius, int number_of_sides, const Vec& circle_normal = Vec(0,0,1), const Vec& circle_center = Vec(0,0,0));
	Circle& operator= (const Circle& from);

	Vec & getNormal();
	Vec & getCenter();

	Vector<Vec> getPoints();

	void translate(const Vec & to);

	void draw(double lineWidth, const Color4 & color);
};

#endif // CIRCLE_H
