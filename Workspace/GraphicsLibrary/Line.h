#ifndef LINE_H
#define LINE_H

#include "Point.h"
#include "Color4.h"

typedef std::pair<double, double> Pairdouble;

class Line
{
public:

	Vec a, b;
	double length;
	int index;

	Color4 color;

	Line();
	Line(const Line&);
	Line(const Vec& from, const Vec& to, int i = -1, 
		const Color4& newColor = Color4(255, 164, 0));
	Line(const Vec& start, const Vec& direction, double length, int i = -1, 
		const Color4& newColor = Color4(255, 164, 0));

	Line colored(const Color4& newColor);

	void set(const Vec& from, const Vec& to);

	Vec direction() const;

	void reverse();

	bool hasPoint(const Vec& point, double eps = Epsilon);

	Vec pointAt(double time) const;
	Vec project(const Vec& point);
	double timeAt(const Vec& point);
	Pairdouble lengthsAt(const Vec& point);
	Pairdouble lengthsAt(double time);
	Vec midPoint();

	void translateBy(const Vec& delta);
        void rotateAroundStart(const qglviewer::Quaternion& q);

	double distanceToUnbounded(const Vec& point);

	void intersectLine( const Line& S2, Vec & pa, Vec & pb );

	void setColor(const Color4& newColor);

	void draw();
};

#endif // LINE_H
