#include "Plane.h"
#include "Utility/SimpleDraw.h"

Vec3d Plane::projectionOf(const Vec3d &point)
{
	if(IsOn(point))
		return point;
	else
	{
		double dist = GetPointDistance(point);
		return point - dist * n;
	}
}

int Plane::LineIntersect(Line l, Vec3d & result)
{
	return LineIntersect(l.a, l.b, center, result);
}

int Plane::LineIntersect(const Vec3d& start, const Vec3d& end, const Vec3d& pointOnPlane, Vec3d & result, double Epsilon)
{
	double denom, numer;

	Vec3d dir = end - start;

	denom = dot(this->n, dir);

	// Check if zero or close to zero
	if (denom > -Epsilon && denom < Epsilon) 
		return PARALLEL;

	numer = dot(this->n, (pointOnPlane - start));

	double u = numer / denom;

	// Is there an intersection along the line?
	if (u < 0.0 || u > 1.0)
	{
		return NO_INTERSECT;
	}
	else
	{
		result = start + (dir * u); // point of intersection

		if(u == 0.0 || u == 1.0)
			return ENDPOINT_INTERSECT;
		else
			return INTERSECT;
	}
}

bool Plane::IsInTri(const Vec3d& p, const Vec3d &a, const Vec3d &b, const Vec3d &c) const
{
	return isSameSide(p, a,b,c) && isSameSide(p, b,a,c) && isSameSide(p, c,a,b);
}

/*
* See http://local.wasp.uwa.edu.au/~pbourke/papers/conrec/ for the original
* paper by Paul D. Bourke.
*/

/*-------------------------------------------------------------------------
Create a contour slice through a 3 vertex facet "p"
Given the normal of the cutting plane "n" and a point on the plane "p0"
Return
0 if the contour plane doesn't cut the facet
2 if it does cut the facet, the contour line segment is p1->p2
-1 for an unexpected occurrence
If a vertex touches the contour plane nothing need to be drawn!?
Note: the following has been written as a "stand along" piece of
code that will work but is far from efficient....
*/

int Plane::ContourFacet( Vec3d a, Vec3d b, Vec3d c, Vec3d & p1, Vec3d & p2 )
{
	double A,B,C,D;
	double side[3];

	Vec3d p [3] = {a, b, c}; 

	/*
	Determine the equation of the plane as
	Ax + By + Cz + D = 0
	*/
	A = n.x();
	B = n.y();
	C = n.z();
	D = -(n.x()*center.x() + n.y()*center.y() + n.z()*center.z());

	/*
	Evaluate the equation of the plane for each vertex
	If side < 0 then it is on the side to be retained
	else it is to be clipped
	*/
	side[0] = A*p[0].x() + B*p[0].y() + C*p[0].z() + D;
	side[1] = A*p[1].x() + B*p[1].y() + C*p[1].z() + D;
	side[2] = A*p[2].x() + B*p[2].y() + C*p[2].z() + D;

	/* Are all the vertices on the same side */
	if (side[0] >= 0 && side[1] >= 0 && side[2] >= 0)	return 0 ;
	if (side[0] <= 0 && side[1] <= 0 && side[2] <= 0)	return 0 ;

	/* Is p0 the only point on a side by itself */
	if ((SIGN(side[0]) != SIGN(side[1])) && (SIGN(side[0]) != SIGN(side[2]))) {
		p1.x() = p[0].x() - side[0] * (p[2].x() - p[0].x()) / (side[2] - side[0]);
		p1.y() = p[0].y() - side[0] * (p[2].y() - p[0].y()) / (side[2] - side[0]);
		p1.z() = p[0].z() - side[0] * (p[2].z() - p[0].z()) / (side[2] - side[0]);
		p2.x() = p[0].x() - side[0] * (p[1].x() - p[0].x()) / (side[1] - side[0]);
		p2.y() = p[0].y() - side[0] * (p[1].y() - p[0].y()) / (side[1] - side[0]);
		p2.z() = p[0].z() - side[0] * (p[1].z() - p[0].z()) / (side[1] - side[0]);
		return(2);
	}

	/* Is p1 the only point on a side by itself */
	if ((SIGN(side[1]) != SIGN(side[0])) && (SIGN(side[1]) != SIGN(side[2]))) {
		p1.x() = p[1].x() - side[1] * (p[2].x() - p[1].x()) / (side[2] - side[1]);
		p1.y() = p[1].y() - side[1] * (p[2].y() - p[1].y()) / (side[2] - side[1]);
		p1.z() = p[1].z() - side[1] * (p[2].z() - p[1].z()) / (side[2] - side[1]);
		p2.x() = p[1].x() - side[1] * (p[0].x() - p[1].x()) / (side[0] - side[1]);
		p2.y() = p[1].y() - side[1] * (p[0].y() - p[1].y()) / (side[0] - side[1]);
		p2.z() = p[1].z() - side[1] * (p[0].z() - p[1].z()) / (side[0] - side[1]);
		return(2);
	}

	/* Is p2 the only point on a side by itself */
	if ((SIGN(side[2]) != SIGN(side[0])) && (SIGN(side[2]) != SIGN(side[1]))) {
		p1.x() = p[2].x() - side[2] * (p[0].x() - p[2].x()) / (side[0] - side[2]);
		p1.y() = p[2].y() - side[2] * (p[0].y() - p[2].y()) / (side[0] - side[2]);
		p1.z() = p[2].z() - side[2] * (p[0].z() - p[2].z()) / (side[0] - side[2]);
		p2.x() = p[2].x() - side[2] * (p[1].x() - p[2].x()) / (side[1] - side[2]);
		p2.y() = p[2].y() - side[2] * (p[1].y() - p[2].y()) / (side[1] - side[2]);
		p2.z() = p[2].z() - side[2] * (p[1].z() - p[2].z()) / (side[1] - side[2]);
		return(2);
	}

	/* Shouldn't get here */
	return 0 ;
}

void Plane::draw(double extent)
{
	glDisable(GL_LIGHTING);

	Vec3d u = orthogonalVector(n).normalized();
	Vec3d v = cross(n, u).normalized();

	u *= extent;
	v *= extent;

	Vec3d p1 = center + (u + v);
	Vec3d p2 = center + (-u + v);
	Vec3d p3 = center + (-v + u);
	Vec3d p4 = center + (-v + -u);

	// Rendering options
	glEnable (GL_POINT_SMOOTH);
	glEnable (GL_LINE_SMOOTH);
	glEnable (GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glHint (GL_LINE_SMOOTH_HINT, GL_NICEST);

	float color[4];
	glGetFloatv(GL_CURRENT_COLOR, color);

	// Draw Borders
	glColor4f(color[0]*0.8, color[1]*0.8, color[2]*0.8, 0.5);
	glLineWidth(2.0);
	glBegin(GL_LINE_STRIP);
	glVertex3dv(p1);
	glVertex3dv(p2);
	glVertex3dv(p4);
	glVertex3dv(p3);
	glVertex3dv(p1);
	glEnd();

	// Draw Center
	glColor3f(color[0], color[1], color[2]);
	glPointSize(4.0);
	glBegin(GL_POINTS);
	glVertex3dv(center);
	glEnd();
	glPointSize(8.0);
	glColor4f(1, 1, 1, 0.5);
	glBegin(GL_POINTS);
	glVertex3dv(center);
	glEnd();

	// Draw Normal
	glBegin(GL_LINES);
	glVertex3dv(center);
	glVertex3dv(center + (n * 0.2 * extent));
	glEnd();

	// Draw Transparent Fills
	glColor4f(color[0], color[1], color[2], 0.05f);
	glBegin(GL_QUADS);
	glVertex3dv(p1);
	glVertex3dv(p2);
	glVertex3dv(p4);
	glVertex3dv(p3);
	glEnd();

	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
}

void Plane::getTangents( Vec3d& X, Vec3d& Y ) const
{
	X = orthogonalVector(n).normalized();
	Y = cross(n, X).normalized();
}

double Plane::GetPointDistance( const Vec3d& pt ) const
{
	return d + dot(n, pt);
}

void Plane::projectLine( Line & line )
{
	Vec3d a = line.a, b = line.b;

	line.a = projectionOf(a);
	line.b = projectionOf(b);

	line.length = (line.a-line.b).norm();
}

Vec3d Plane::reflection( const Vec3d& v )
{
	return n * GetPointDistance(v) * -2.0 + v;
}
