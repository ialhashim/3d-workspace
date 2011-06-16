#include "Circle.h"

Circle::Circle()
{
	this->radius = 1.0;
	this->numSides = 0;
}

Circle::Circle(float from_radius, int number_of_sides, const Vec& circle_normal, const Vec& circle_center)
{
	this->radius = from_radius;
	this->numSides = number_of_sides;
	this->normal = circle_normal;

	const Vec axis = normal;
	Vec ref = axis.orthogonalVec();

	double theta = (2.0 * M_PI) / (double) numSides;
	double startAngle = 0;

	for(int i = 0; i < numSides; i++)
	{
                point.push_back(radius * qglviewer::Quaternion(axis, startAngle + theta).rotate(ref).unit());
		startAngle += theta;
	}

	translate(circle_center);
}

Circle& Circle::operator= (const Circle& from)
{
	this->radius = from.radius;
	this->normal = from.normal;
	this->numSides = from.numSides;
	this->point = from.point;
	this->center = from.center;

	return *this;
}

void Circle::translate(const Vec & to)
{
	Vec delta = to - center;

	for(int i = 0; i < (int)point.size(); i++)
		point[i] += delta;

	center = to;
}

Vec & Circle::getNormal()
{
	return normal;
}

Vec & Circle::getCenter()
{
	return center;
}

Vector<Vec> Circle::getPoints()
{
	return point;
}

void Circle::draw(double lineWidth, const Color4 & color)
{
	glLineWidth(lineWidth);
	glColor4ubv(color.m_v);

	glBegin(GL_LINE_STRIP);
	for(unsigned int i = 0; i <= point.size(); i++)
	{
		glVertex3fv(point[i % point.size()]);

	}
	glEnd();
}
