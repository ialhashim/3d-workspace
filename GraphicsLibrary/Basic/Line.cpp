#include "Line.h"
#include "Utility/SimpleDraw.h"

Line::Line()
{
	length = 0;
	index = -1;
}

Line::Line(const Line& from)
{
	this->a = from.a;
	this->b = from.b;
	this->length = from.length;
	this->index = from.index;
	this->color = from.color;
}

Line::Line( const Vec3d& from, const Vec3d& to, int i, const Color& newColor)
{
	a = from;
	b = to;

	length = (a-b).norm();

	index = i;

	color = newColor;
}

Line::Line( const Vec3d& start, const Vec3d& direction, double length, int i, const Color& newColor)
{
	a = start;
	b = start + (direction.normalized() * length);

	length = (a-b).norm();

	index = i;

	color = newColor;
}

void Line::set( const Vec3d& from, const Vec3d& to )
{
	a = from;
	b = to;
	length = (a-b).norm();
	index = -1;
	color = Color(1, 0.5, 0, 1);
}

Line Line::colored( const Color& newColor )
{
	Line l (*this);
	l.setColor(newColor);
	return l;
}

void Line::reverse()
{
	Vec3d temp = a;
	a = b;
	b = temp;
}

Vec3d Line::direction() const
{
	return b - a;
}

bool Line::hasPoint( const Vec3d& point, double eps)
{
	double dist = (point - a).norm() + (point - b).norm();

	if(dist < length + eps)
		return true;
	else
		return false;
}

double Line::distanceToUnbounded( const Vec3d& point )
{
	Vec3d O = a, P = point, D = direction().normalized();

	Vec3d closest = O + (dot((P-O),D))/dot(D,D) * D;

	return (closest - point).norm();
}

Vec3d Line::midPoint()
{
	return (a + b) * 0.5;
}

Vec3d Line::project( const Vec3d& point )
{
	return pointAt(dot((point - a) , direction()));
}

Vec3d Line::pointAt( double time ) const
{
	double dist = time * length;
	return a + (direction().normalized() * dist);
}

double Line::timeAt( const Vec3d& point )
{
	return (point - a).norm() / length;
}

Pairdouble Line::lengthsAt( const Vec3d& point )
{
	double dist1 = (point - a).norm();
	double dist2 = length - dist1;

	return Pairdouble(dist1, dist2);
}

Pairdouble Line::lengthsAt( double time )
{
	time = Max(0, Min(1.0, time)); // bounded
	return Pairdouble(length * time, length * (1.0 - time));
}

void Line::translateBy( const Vec3d& delta )
{
	a += delta;
	b += delta;
}

void Line::intersectLine( const Line& S2, Vec3d & pa, Vec3d & pb, double Epsilon )
{
	double EPS = Epsilon; // experimental

	Vec3d   u = this->b - this->a;
	Vec3d   v = S2.b - S2.a;
	Vec3d   w = this->a - S2.a;

	double    k = dot(u,u);			// (a) always >= 0
	double    j = dot(u,v);			// (b)
	double    c = dot(v,v);			// always >= 0
	double    d = dot(u,w);
	double    e = dot(v,w);
	double    D = k*c - j*j;       // always >= 0
	double    sc, sN, sD = D;      // sc = sN / sD, default sD = D >= 0
	double    tc, tN, tD = D;      // tc = tN / tD, default tD = D >= 0

	// compute the line parameters of the two closest points
	if (D < EPS) { // the lines are almost parallel
		sN = 0.0;        // force using point P0 on segment S1
		sD = 1.0;        // to prevent possible division by 0.0 later
		tN = e;
		tD = c;
	}
	else {                // get the closest points on the infinite lines
		sN = (j*e - c*d);
		tN = (k*e - j*d);
		if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
			sN = 0.0;
			tN = e;
			tD = c;
		}
		else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
			sN = sD;
			tN = e + j;
			tD = c;
		}
	}

	if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
		tN = 0.0;
		// recompute sc for this edge
		if (-d < 0.0)
			sN = 0.0;
		else if (-d > k)
			sN = sD;
		else {
			sN = -d;
			sD = k;
		}
	}
	else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
		tN = tD;
		// recompute sc for this edge
		if ((-d + j) < 0.0)
			sN = 0;
		else if ((-d + j) > k)
			sN = sD;
		else {
			sN = (-d + j);
			sD = k;
		}
	}
	// finally do the division to get sc and tc
	sc = (abs(sN) < EPS ? 0.0 : sN / sD);
	tc = (abs(tN) < EPS ? 0.0 : tN / tD);

	pa = this->pointAt(sc);
	pb = S2.pointAt(tc);
}

void Line::draw()
{
	// Rendering options
	glDisable (GL_LIGHTING);
	glEnable (GL_POINT_SMOOTH);
	glEnable (GL_LINE_SMOOTH);
	glEnable (GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glHint (GL_LINE_SMOOTH_HINT, GL_NICEST);

	glColor3f(color[0], color[1], color[2]);

	glLineWidth(2.0f);

	glBegin(GL_LINES);
	glVertex3dv(a);
	glVertex3dv(b);
	glEnd();

	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
}

void Line::setColor( const Color& newColor )
{
	color = newColor;
}

std::vector<Vec3d> Line::uniformSample( int numSamples )
{
	std::vector<Vec3d> result;

	double deltaLength = length / (numSamples - 1);
	Vec3d delta = deltaLength * direction().normalized();

	result.push_back(a);

	for(uint i = 1; i < numSamples; i++)
		result.push_back(result.back() + delta);

	return result;
}

Line::operator const std::vector<Vec3d>()
{
	return uniformSample(2);
}
