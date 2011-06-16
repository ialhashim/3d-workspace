#include "Plane.h"

Vec Plane::projectionOf(const Vec &point)
{
	if(IsOn(point))
		return point;
	else
	{
		Vec proj = point - this->center;

		proj.projectOnPlane(n);

		return proj + this->center;
	}
}

int Plane::LineIntersect(const Line& l, Vec & result)
{
    return LineIntersect(l.a, l.b, center, result);
}

int Plane::LineIntersect(const Vec& start, const Vec& end, const Vec& pointOnPlane, Vec & result)
{
	double denom, numer;

	Vec dir = end - start;

	denom = this->n * dir;

	// Check if zero or close to zero
	if (denom > -Epsilon && denom < Epsilon) 
		return PARALLEL;

	numer = this->n * (pointOnPlane - start);

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

bool Plane::IsInTri(const Vec& p, const Vec &a, const Vec &b, const Vec &c) const
{
	return Point3D::isSameSide(p,a, b,c) 
		&& Point3D::isSameSide(p,b, a,c) 
		&& Point3D::isSameSide(p,c, a,b);
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
int Plane::ContourFacet( BaseTriangle * f, Vec & p1, Vec & p2 )
{
   double A,B,C,D;
   double side[3];

   Vec p [3] = {f->vec(0), f->vec(1), f->vec(2)}; 

   /*
      Determine the equation of the plane as
      Ax + By + Cz + D = 0
   */
   A = n.x;
   B = n.y;
   C = n.z;
   D = -(n.x*center.x + n.y*center.y + n.z*center.z);

   /*
      Evaluate the equation of the plane for each vertex
      If side < 0 then it is on the side to be retained
      else it is to be clipped
   */
   side[0] = A*p[0].x + B*p[0].y + C*p[0].z + D;
   side[1] = A*p[1].x + B*p[1].y + C*p[1].z + D;
   side[2] = A*p[2].x + B*p[2].y + C*p[2].z + D;

   /* Are all the vertices on the same side */
   if (side[0] >= 0 && side[1] >= 0 && side[2] >= 0)	return 0 ;
   if (side[0] <= 0 && side[1] <= 0 && side[2] <= 0)	return 0 ;

   /* Is p0 the only point on a side by itself */
   if ((SIGN(side[0]) != SIGN(side[1])) && (SIGN(side[0]) != SIGN(side[2]))) {
      p1.x = p[0].x - side[0] * (p[2].x - p[0].x) / (side[2] - side[0]);
      p1.y = p[0].y - side[0] * (p[2].y - p[0].y) / (side[2] - side[0]);
      p1.z = p[0].z - side[0] * (p[2].z - p[0].z) / (side[2] - side[0]);
      p2.x = p[0].x - side[0] * (p[1].x - p[0].x) / (side[1] - side[0]);
      p2.y = p[0].y - side[0] * (p[1].y - p[0].y) / (side[1] - side[0]);
      p2.z = p[0].z - side[0] * (p[1].z - p[0].z) / (side[1] - side[0]);
      return(2);
   }

   /* Is p1 the only point on a side by itself */
   if ((SIGN(side[1]) != SIGN(side[0])) && (SIGN(side[1]) != SIGN(side[2]))) {
      p1.x = p[1].x - side[1] * (p[2].x - p[1].x) / (side[2] - side[1]);
      p1.y = p[1].y - side[1] * (p[2].y - p[1].y) / (side[2] - side[1]);
      p1.z = p[1].z - side[1] * (p[2].z - p[1].z) / (side[2] - side[1]);
      p2.x = p[1].x - side[1] * (p[0].x - p[1].x) / (side[0] - side[1]);
      p2.y = p[1].y - side[1] * (p[0].y - p[1].y) / (side[0] - side[1]);
      p2.z = p[1].z - side[1] * (p[0].z - p[1].z) / (side[0] - side[1]);
      return(2);
   }

   /* Is p2 the only point on a side by itself */
   if ((SIGN(side[2]) != SIGN(side[0])) && (SIGN(side[2]) != SIGN(side[1]))) {
      p1.x = p[2].x - side[2] * (p[0].x - p[2].x) / (side[0] - side[2]);
      p1.y = p[2].y - side[2] * (p[0].y - p[2].y) / (side[0] - side[2]);
      p1.z = p[2].z - side[2] * (p[0].z - p[2].z) / (side[0] - side[2]);
      p2.x = p[2].x - side[2] * (p[1].x - p[2].x) / (side[1] - side[2]);
      p2.y = p[2].y - side[2] * (p[1].y - p[2].y) / (side[1] - side[2]);
      p2.z = p[2].z - side[2] * (p[1].z - p[2].z) / (side[1] - side[2]);
      return(2);
   }

   /* Shouldn't get here */
   return 0 ;
}

void Plane::draw(double extent)
{
	glDisable(GL_LIGHTING);

	Vec u = (n.orthogonalVec()).unit();
	Vec v = (n ^ u).unit();

	u *= extent;
	v *= extent;

	Vec p1 = center + (u + v);
	Vec p2 = center + (-u + v);
	Vec p3 = center + (-v + u);
	Vec p4 = center + (-v + -u);

	// Rendering options
	glEnable (GL_POINT_SMOOTH);
	glEnable (GL_LINE_SMOOTH);
	glEnable (GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glHint (GL_LINE_SMOOTH_HINT, GL_NICEST);

	float color[4];
	glGetFloatv(GL_CURRENT_COLOR, color);

	// Draw Borders
	glColor3f(color[0]*0.8, color[1]*0.8, color[2]*0.8);
	glLineWidth(3.0);
	glBegin(GL_LINE_STRIP);
		glVertex3fv(p1);
		glVertex3fv(p2);
		glVertex3fv(p4);
		glVertex3fv(p3);
		glVertex3fv(p1);
	glEnd();

	// Draw Center
	glColor3f(color[0], color[1], color[2]);
	glPointSize(8.0);
	glBegin(GL_POINTS);
		glVertex3fv(center);
	glEnd();
	glPointSize(12.0);
	glColor3f(1, 1, 1);
	glBegin(GL_POINTS);
		glVertex3fv(center);
	glEnd();

	// Draw Normal
	glBegin(GL_LINES);
		glVertex3fv(center);
		glVertex3fv(center + n);
	glEnd();

	// Draw Transpernt Fills
	glColor4f(color[0], color[1], color[2], 0.15f);
	glBegin(GL_QUADS);
		glVertex3fv(p1);
		glVertex3fv(p2);
		glVertex3fv(p4);
		glVertex3fv(p3);
	glEnd();

	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
}

void Plane::getTangents( Vec& X, Vec& Y ) const
{
	X = n.orthogonalVec().unit();
	Y = (n ^ X).unit();
}

double Plane::GetPointDistance( const Vec& pt ) const
{
	return d + n * pt;
}

void Plane::projectLine( Line & line )
{
        Vec a = line.a, b = line.b;

        line.a = projectionOf(a);
        line.b = projectionOf(b);

        line.length = (line.a-line.b).norm();
}
