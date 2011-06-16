#include "Line.h"

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

Line::Line( const Vec& from, const Vec& to, int i, const Color4& newColor)
{
	a = from;
	b = to;

	length = (a-b).norm();

	index = i;

	color = newColor;
}

Line::Line( const Vec& start, const Vec& direction, double length, int i, const Color4& newColor)
{
	a = start;
	b = start + (direction.unit() * length);

	length = (a-b).norm();

	index = i;

	color = newColor;
}

void Line::set( const Vec& from, const Vec& to )
{
	a = from;
	b = to;
	length = (a-b).norm();
	index = -1;
	color = Color4(255, 164, 0);
}

Line Line::colored( const Color4& newColor )
{
	Line l (*this);
	l.setColor(newColor);
	return l;
}

void Line::reverse()
{
	Vec temp = a;
	a = b;
	b = temp;
}

Vec Line::direction() const
{
	return b - a;
}

bool Line::hasPoint( const Vec& point, double eps)
{
	double dist = (point - a).norm() + (point - b).norm();

	if(dist < length + eps)
		return true;
	else
		return false;
}

double Line::distanceToUnbounded( const Vec& point )
{
	Vec O = a, P = point, D = direction().unit();

	Vec closest = O + ((P-O)*D)/(D*D) * D;

	return (closest - point).norm();
}

Vec Line::midPoint()
{
	return (a + b) * 0.5f;
}

Vec Line::project( const Vec& point )
{
	return pointAt((point - a) * direction());
}

Vec Line::pointAt( double time ) const
{
	double dist = time * length;
	return a + (direction().unit() * dist);
}

double Line::timeAt( const Vec& point )
{
	return (point - a).norm() / length;
}

Pairdouble Line::lengthsAt( const Vec& point )
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

void Line::translateBy( const Vec& delta )
{
	a += delta;
	b += delta;
}

void Line::rotateAroundStart( const qglviewer::Quaternion& q )
{
	b = Point3D::RotateAround(b, a, q);
}

/*bool Line::intersectLine( Line * other, Vec * pa, Vec * pb )
{
Vec p13, p43, p21;

double d1343,d4321,d1321,d4343,d2121;
double numer,denom;

p13 = this->a - other->a;
p43 = other->b - other->a;

double EPS = 1.0e-12f;

if (abs(p43.x)  < EPS && abs(p43.y)  < EPS && abs(p43.z)  < EPS)
return(FALSE);

p21 = this->b - this->a;

if (abs(p21.x)  < EPS && abs(p21.y)  < EPS && abs(p21.z)  < EPS)
return(FALSE);

d1343 = p13.x * p43.x + p13.y * p43.y + p13.z * p43.z;
d4321 = p43.x * p21.x + p43.y * p21.y + p43.z * p21.z;
d1321 = p13.x * p21.x + p13.y * p21.y + p13.z * p21.z;
d4343 = p43.x * p43.x + p43.y * p43.y + p43.z * p43.z;
d2121 = p21.x * p21.x + p21.y * p21.y + p21.z * p21.z;

denom = d2121 * d4343 - d4321 * d4321;

if (abs(denom) < EPS)
return(FALSE);

numer = d1343 * d4321 - d1321 * d4343;

double mua = numer / denom;
double mub = (d1343 + d4321 * (mua)) / d4343;

if(mua < 0 || mua > 1 || mub < 0 || mub > 1)
return false;

pa->x = a.x + mua * p21.x;
pa->y = a.y + mua * p21.y;
pa->z = a.z + mua * p21.z;

pb->x = other->a.x + mub * p43.x;
pb->y = other->a.y + mub * p43.y;
pb->z = other->a.z + mub * p43.z;

return(TRUE);
}*/

void Line::intersectLine( const Line& S2, Vec & pa, Vec & pb )
{
	double EPS = Epsilon; // experimental

	Vec   u = this->b - this->a;
	Vec   v = S2.b - S2.a;
	Vec   w = this->a - S2.a;

	float    k = (u*u);			// (a) always >= 0
	float    j = (u*v);			// (b)
	float    c = (v*v);			// always >= 0
	float    d = (u*w);
	float    e = (v*w);
	float    D = k*c - j*j;       // always >= 0
	float    sc, sN, sD = D;      // sc = sN / sD, default sD = D >= 0
	float    tc, tN, tD = D;      // tc = tN / tD, default tD = D >= 0

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

	glColor3f(color.r(), color.g(), color.b());

	glLineWidth(2.0f);

	glBegin(GL_LINES);
	glVertex3fv(a);
	glVertex3fv(b);
	glEnd();

	/*glPointSize(12.0f);
	glBegin(GL_POINTS);
	glVertex3fv(b);
	glColor3f(color.r() *0.5f, color.g() *0.5f, color.b() *0.5f); // darker
	glVertex3fv(a);
	glEnd();*/

	/*glColor3f(1,1,1);
	glPointSize(14.0f);
	glBegin(GL_POINTS);
	glVertex3fv(a);
	glVertex3fv(b);
	glEnd();*/

	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
}

void Line::setColor( const Color4& newColor )
{
	color = newColor;
}
