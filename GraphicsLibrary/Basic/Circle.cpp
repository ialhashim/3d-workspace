#include "Circle.h"

#include "GUI/Viewer/libQGLViewer/QGLViewer/qglviewer.h"

Circle::Circle()
{
	this->radius = 1.0;
	this->numSides = 0;
}

Circle::Circle( const Vec3d& circle_center, const Vec3d& circle_normal, double from_radius /*= 1.0*/, int number_of_sides /*= 2.0*/ )
{
	this->radius = from_radius;
	this->numSides = number_of_sides;
	this->normal = circle_normal;
	this->center = circle_center;

	Vec3d refV(V2V(qglviewer::Vec(normal[0],normal[1],normal[2]).orthogonalVec()));

	double theta = (2.0 * M_PI) / (double) numSides;
	double startAngle = M_PI / 4;

	for(int i = 0; i < numSides; i++)
	{
		Vec3d v = ROTATED_VEC(refV, startAngle + theta, normal);
        point.push_back(center + (v * radius));

		startAngle += theta;
	}
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

void Circle::translate(const Vec3d & to)
{
	Vec3d delta = to - center;

	for(int i = 0; i < (int)point.size(); i++)
		point[i] += delta;

	center = to;
}

Vec3d & Circle::getNormal()
{
	return normal;
}

Vec3d & Circle::getCenter()
{
	return center;
}

StdVector<Vec3d> Circle::getPoints()
{
	return point;
}

void Circle::draw(double lineWidth, const Vec4d & color)
{
	glLineWidth(lineWidth);
	glColor4dv(color);

	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable (GL_LINE_SMOOTH);

	glBegin(GL_LINE_LOOP);
	for(int i = 0; i <= point.size(); i++)
		glVertex3dv(point[i % point.size()]);
	glEnd();

	// Fill empty spaces
	//glPointSize(lineWidth);
	//glBegin(GL_POINTS);
	//for(int i = 0; i < point.size(); i++) glVertex3dv(point[i]);
	//glEnd();
}

void Circle::drawFilled(const Vec4d & fillColor, double lineWidth, const Vec4d & borderColor)
{
	draw(lineWidth, borderColor);

	glColor4dv(fillColor);

	glBegin(GL_TRIANGLE_FAN);
	glVertex3dv(center);
	for(int i = 0; i <= point.size(); i++)
	{
		glVertex3dv(point[i % point.size()]);
	}
	glEnd();
}
