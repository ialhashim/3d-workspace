#pragma once

#include "Spline.h"

#include "BezierSpline.h"

#include "Plane.h"

class UserCurve
{
private:

public:
	UserCurve();

	Spline spline;

	Vec plane_origin;
	Vec plane_normal;

	Vector<Vec> points;

	bool isReady;
	bool isSimplified;
	bool isVisible;

	void setPlane(const Vec & pointOnPlane, const Vec & normal);
	Plane getPlane();
	void addPoint(const Vec & origin, const Vec & direction);
	void moveHead(const Vec & to);

	Vector<Vec> getPath(float stepLength);
	Spline * getSpline();

	static Vec intersectionRayPlane(const Vec & planeNormal, 
		const Vec & planeOrigin, const Vec & rayOrigin, const Vec & rayDirection);

	void simplify();
	void simplify2();

	void draw();
};
