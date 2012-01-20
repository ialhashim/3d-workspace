#ifndef LINE_H
#define LINE_H

#include "Surface_mesh.h"
#include "Utility/Macros.h"

typedef std::pair<double, double> Pairdouble;

class Line
{
public:

	Vec3d a, b;
	double length;
	int index;

	Color color;

	Line();
	Line(const Line&);
	Line(const Vec3d& from, const Vec3d& to, int i = -1, const Color& newColor = Color(1, 0.5, 0, 1));
	Line(const Vec3d& start, const Vec3d& direction, double length, int i = -1, const Color& newColor = Color(1, 0.5, 0, 1));

	Line colored(const Color& newColor);

	void set(const Vec3d& from, const Vec3d& to);

	Vec3d direction() const;

	void reverse();

	bool hasPoint(const Vec3d& point, double eps = 1e-10);

	Vec3d pointAt(double time) const;
	Vec3d project(const Vec3d& point);
	double timeAt(const Vec3d& point);
	Pairdouble lengthsAt(const Vec3d& point);
	Pairdouble lengthsAt(double time);
	Vec3d midPoint();

	std::vector<Vec3d> uniformSample(int numSamples);
	
	operator const std::vector<Vec3d>();

	void translateBy(const Vec3d& delta);

	double distanceToUnbounded(const Vec3d& point);

	void intersectLine( const Line& S2, Vec3d & pa, Vec3d & pb, double Epsilon = 1e-10 );
	void setColor(const Color& newColor);

	void draw();
};

#endif // LINE_H
