#include "SimpleDraw.h"
#include "Macros.h"
#include <GL/glu.h>

// Bad includes.. needed for rotations for now
#include "GUI/Viewer/libQGLViewer/QGLViewer/qglviewer.h"
#include "ColorMap.h"

void SimpleDraw::DrawBox(const Vec3d& center, float width, float length, float height, float r, float g, float b, float lineWidth)
{
	glDisable(GL_LIGHTING);

	glColor3f(r, g, b);

	Vec3d  c1, c2, c3, c4;
	Vec3d  bc1, bc2, bc3, bc4;

	c1 = Vec3d (width, length, height) + center;
	c2 = Vec3d (-width, length, height) + center;
	c3 = Vec3d (-width, -length, height) + center;
	c4 = Vec3d (width, -length, height) + center;

	bc1 = Vec3d (width, length, -height) + center;
	bc2 = Vec3d (-width, length, -height) + center;
	bc3 = Vec3d (-width, -length, -height) + center;
	bc4 = Vec3d (width, -length, -height) + center;

	glLineWidth(lineWidth);

	glBegin(GL_LINES);
	glVertex3dv(c1);glVertex3dv(bc1);
	glVertex3dv(c2);glVertex3dv(bc2);
	glVertex3dv(c3);glVertex3dv(bc3);
	glVertex3dv(c4);glVertex3dv(bc4);
	glVertex3dv(c1);glVertex3dv(c2);
	glVertex3dv(c3);glVertex3dv(c4);
	glVertex3dv(c1);glVertex3dv(c4);
	glVertex3dv(c2);glVertex3dv(c3);
	glVertex3dv(bc1);glVertex3dv(bc2);
	glVertex3dv(bc3);glVertex3dv(bc4);
	glVertex3dv(bc1);glVertex3dv(bc4);
	glVertex3dv(bc2);glVertex3dv(bc3);
	glEnd();

	glEnable(GL_LIGHTING);

}

void SimpleDraw::DrawSolidBox(const Vec3d & center, float width, float length, float height, float r, float g, float b, float a)
{
	glColor3f(r, g, b);

	Vec3d  c1, c2, c3, c4;
	Vec3d  bc1, bc2, bc3, bc4;

	width *= 0.5;
	length *= 0.5;
	height *= 0.5;

	c1 = Vec3d (width, length, height) + center;
	c2 = Vec3d (-width, length, height) + center;
	c3 = Vec3d (-width, -length, height) + center;
	c4 = Vec3d (width, -length, height) + center;

	bc1 = Vec3d (width, length, -height) + center;
	bc2 = Vec3d (-width, length, -height) + center;
	bc3 = Vec3d (-width, -length, -height) + center;
	bc4 = Vec3d (width, -length, -height) + center;

	glShadeModel(GL_FLAT);

	SimpleDraw::DrawSquare(c1, c2, c3, c4, true, r,g,b,a);
	SimpleDraw::DrawSquare(bc4, bc3, bc2, bc1, true, r,g,b,a);

	SimpleDraw::DrawSquare(c1, c4, bc4, bc1, true, r,g,b,a);
	SimpleDraw::DrawSquare(c2, c1, bc1, bc2, true, r,g,b,a);

	SimpleDraw::DrawSquare(c4, c3, bc3, bc4, true, r,g,b,a);
	SimpleDraw::DrawSquare(c2, bc2, bc3, c3, true, r,g,b,a);

	glShadeModel(GL_SMOOTH);
}

void SimpleDraw::DrawTriangle(const Vec3d & v1, const Vec3d & v2, const Vec3d & v3, float r, float g, float b, float a, bool isOpaque)
{
	glDisable(GL_LIGHTING);

	if(a < 1.0f)
	{
		glEnable(GL_BLEND); 
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}

	// Draw the filled triangle
	if(isOpaque)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset( 0.5, 0.5 );

		glColor4f(r * 0.75f, g * 0.75f, b * 0.75f, RANGED(0,a,1.0f));

		if(a < 1.0f){
			glEnable(GL_BLEND); 
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		}

		glBegin(GL_TRIANGLES);
		glVertex3dv(v1);
		glVertex3dv(v2);
		glVertex3dv(v3);
		glEnd();

		glDisable(GL_POLYGON_OFFSET_FILL);
	}

	// Draw the edges
	glLineWidth(1.0f);

	glColor3f(r * 0.5f, g * 0.5f, b * 0.5f);

	glBegin(GL_LINE_STRIP);
	glVertex3dv(v1);
	glVertex3dv(v2);
	glVertex3dv(v3);
	glVertex3dv(v1);
	glEnd();

	// Draw the points
	int pointSize = 4;

	glEnable(GL_POINT_SMOOTH);

	// Colored dot
	glColor3f(r,g,b);
	glPointSize(pointSize);
	glBegin(GL_POINTS);
	glVertex3dv(v1);
	glVertex3dv(v2);
	glVertex3dv(v3);
	glEnd();

	// White Border
	glPointSize(pointSize + 2);
	glColor3f(1, 1, 1);

	glBegin(GL_POINTS);
	glVertex3dv(v1);
	glVertex3dv(v2);
	glVertex3dv(v3);
	glEnd();

	glDisable (GL_BLEND);

	glEnable(GL_LIGHTING);
}

void SimpleDraw::DrawTriangles( const StdVector< StdVector<Vec3d> > & tris, float r, float g, float b, float a, bool isOpaque, bool isDrawVec3ds)
{
	glDisable(GL_LIGHTING);

	if(a < 1.0f)
	{
		glEnable(GL_BLEND); 
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}

	// Draw the filled triangle
	if(isOpaque)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset( 0.5, 0.5 );

		glColor4f(r * 0.75f, g * 0.75f, b * 0.75f, RANGED(0,a,1.0f));

		if(a < 1.0f){
			glEnable(GL_BLEND); 
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		}

		glBegin(GL_TRIANGLES);
		for(int i = 0; i < (int)tris.size(); i++)
		{
			glVertex3dv(tris[i][0]);
			glVertex3dv(tris[i][1]);
			glVertex3dv(tris[i][2]);
		}
		glEnd();

		glDisable(GL_POLYGON_OFFSET_FILL);
	}

	// Draw the edges
	glLineWidth(1.0f);

	if(isOpaque)
		glColor3f(r * 0.5f, g * 0.5f, b * 0.5f);
	else
		glColor3f(r, g, b);

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	glBegin(GL_TRIANGLES);
	for(int i = 0; i < (int)tris.size(); i++)
	{
		glVertex3dv(tris[i][0]);
		glVertex3dv(tris[i][1]);
		glVertex3dv(tris[i][2]);
	}
	glEnd();

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	if(isDrawVec3ds)
	{
		// Draw the points
		int pointSize = 4;

		glEnable(GL_POINT_SMOOTH);

		// Colored dot
		glColor3f(r,g,b);
		glPointSize(pointSize);
		glBegin(GL_POINTS);
		for(int i = 0; i < (int)tris.size(); i++)
		{
			glVertex3dv(tris[i][0]);
			glVertex3dv(tris[i][1]);
			glVertex3dv(tris[i][2]);
		}
		glEnd();

		// White Border
		glPointSize(pointSize + 2);
		glColor3f(1, 1, 1);

		glBegin(GL_POINTS);
		for(int i = 0; i < (int)tris.size(); i++)
		{
			glVertex3dv(tris[i][0]);
			glVertex3dv(tris[i][1]);
			glVertex3dv(tris[i][2]);
		}
		glEnd();
	}

	glDisable (GL_BLEND);

	glEnable(GL_LIGHTING);
}

void SimpleDraw::DrawPoly( const std::vector<Vec3d> & poly, float r, float g, float b)
{	
	glDisable(GL_LIGHTING);

	glColor3f(r, g, b);

	float lineWidth = 15.0f;
	glLineWidth(lineWidth);

	glBegin(GL_LINE_STRIP);
	for(uint i = 0; i <= poly.size(); i++)
		glVertex3dv(poly[i%poly.size()]);
	glEnd();

	glEnable(GL_LIGHTING);
}

void SimpleDraw::DrawLineTick(const StdVector<Vec3d >& start, const StdVector<Vec3d >& direction, 
							  float len, bool border, float r, float g, float b, float a)
{
	glPushAttrib( GL_ALL_ATTRIB_BITS );

	glDisable(GL_LIGHTING);

	if(a < 1.0f){
		glEnable(GL_BLEND); 
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}

	if(border)
	{
		glClearStencil(0);
		glClear( GL_STENCIL_BUFFER_BIT );
		glEnable( GL_STENCIL_TEST );

		glStencilFunc( GL_ALWAYS, 1, 0xFFFF );
		glStencilOp( GL_KEEP, GL_KEEP, GL_REPLACE );

		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
		glColor3f( 0.0f, 0.0f, 0.0f );
	}

	glLineWidth(1.0f);
	glBegin(GL_LINES);
	glColor4f(r, g, b, a);
	for(int i = 0; i < (int)start.size(); i++)
	{
		glVertex3dv(start[i]);
		glVertex3dv(start[i] + (direction[i] * len));
	}
	glEnd();

	if(border)
	{
		glStencilFunc( GL_NOTEQUAL, 1, 0xFFFF );
		glStencilOp( GL_KEEP, GL_KEEP, GL_REPLACE );

		glLineWidth(4.0f);
		glBegin(GL_LINES);
		glColor4f(0, 0, 0, 1);
		for(int i = 0; i < (int)start.size(); i++)
		{
			glVertex3dv(start[i]);
			glVertex3dv(start[i] + (direction[i] * len * 1.05f));
		}
		glEnd();
	}


	glPopAttrib();

}

void SimpleDraw::DrawSquare(const std::vector<Vec3d> & v, bool isOpaque, float lineWidth, Vec4d color)
{
	DrawSquare(v[0],v[1],v[2],v[3], isOpaque, lineWidth, color[0], color[1], color[2], color[3]);
}

void SimpleDraw::DrawSquare(const Vec3d & v1, const Vec3d & v2, const Vec3d & v3, const Vec3d & v4, 
	bool isOpaque, float lineWidth, Vec4d color)
{
	DrawSquare(v1,v2,v3,v4, isOpaque, lineWidth, color[0], color[1], color[2], color[3]);
}

void SimpleDraw::DrawSquare(const Vec3d & v1, const Vec3d & v2, const Vec3d & v3, const Vec3d & v4, 
							bool isOpaque, float lineWidth, float r, float g, float b, float a)
{
	glEnable(GL_LIGHTING);

	if(isOpaque)
	{
		// Draw the filled square
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset( 0.5f, 0.5f );

		glColor4f(r, g, b, RANGED(0, a, 1.0f));

		glEnable(GL_BLEND); 
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		glBegin(GL_QUADS);
		
		Vec3d v21 = v2 - v1;
		Vec3d v31 = v3 - v1;

		glNormal3dv(cross(v21 , v31).normalized());
		glVertex3dv(v1);
		glVertex3dv(v2);
		glVertex3dv(v3);
		glVertex3dv(v4);
		glEnd();

		glDisable(GL_POLYGON_OFFSET_FILL);
	}

	glDisable(GL_LIGHTING);

	// Draw the edges
	glLineWidth(lineWidth);

	glColor4f(r, g, b, a);

	glBegin(GL_LINE_STRIP);
	glVertex3dv(v1);
	glVertex3dv(v2);
	glVertex3dv(v3);
	glVertex3dv(v4);
	glVertex3dv(v1);
	glEnd();

	glEnable(GL_LIGHTING);
}

void SimpleDraw::DrawSquares( const StdVector<StdVector<Vec3d > >& squares, bool isOpaque, float r, float g, float b, float a )
{
	glDisable(GL_LIGHTING);

	if(isOpaque)
	{
		// Draw the filled square
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset( 0.5, 0.5 );

		glColor4f(r * 0.75f, g * 0.75f, b * 0.75f, RANGED(0, a, 1.0f));

		if(a < 1.0f){
			glEnable(GL_BLEND); 
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		}

		glBegin(GL_QUADS);
		for(int i = 0; i < (int)squares.size(); i++)
		{
			glVertex3dv(squares[i][0]);
			glVertex3dv(squares[i][1]);
			glVertex3dv(squares[i][2]);
			glVertex3dv(squares[i][3]);
		}
		glEnd();

		glDisable(GL_POLYGON_OFFSET_FILL);
	}

	// Draw the edges
	glLineWidth(1.0f);

	glColor3f(r, g, b);

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	glBegin(GL_QUADS);
	for(int i = 0; i < (int)squares.size(); i++)
	{
		glVertex3dv(squares[i][0]);
		glVertex3dv(squares[i][1]);
		glVertex3dv(squares[i][2]);
		glVertex3dv(squares[i][3]);
	}
	glEnd();

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glEnable(GL_LIGHTING);
}

void SimpleDraw::DrawCube( const Vec3d & center, float length /*= 1.0f*/ )
{
	static GLdouble n[6][3] ={
		{-1.0, 0.0, 0.0},
		{0.0, 1.0, 0.0},
		{1.0, 0.0, 0.0},
		{0.0, -1.0, 0.0},
		{0.0, 0.0, 1.0},
		{0.0, 0.0, -1.0}};

	static GLint faces[6][4] ={
		{0, 1, 2, 3},
		{3, 2, 6, 7},
		{7, 6, 5, 4},
		{4, 5, 1, 0},
		{5, 6, 2, 1},
		{7, 4, 0, 3}};

	GLdouble v[8][3];GLint i;

	v[0][0] = v[1][0] = v[2][0] = v[3][0] = -length / 2;
	v[4][0] = v[5][0] = v[6][0] = v[7][0] = length / 2;
	v[0][1] = v[1][1] = v[4][1] = v[5][1] = -length / 2;
	v[2][1] = v[3][1] = v[6][1] = v[7][1] = length / 2;
	v[0][2] = v[3][2] = v[4][2] = v[7][2] = -length / 2;
	v[1][2] = v[2][2] = v[5][2] = v[6][2] = length / 2;

	glPushMatrix();
	glTranslatef(center.x(), center.y(), center.z());

	for (i = 0; i < 6; i++) 
	{
		glBegin(GL_QUADS);
		glNormal3dv(&n[i][0]);
		glVertex3dv(&v[faces[i][0]][0]);
		glVertex3dv(&v[faces[i][1]][0]);
		glVertex3dv(&v[faces[i][2]][0]);
		glVertex3dv(&v[faces[i][3]][0]);
		glEnd();
	}

	glPopMatrix();

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void SimpleDraw::DrawSphere( const Vec3d & center, float radius /*= 1.0f*/ )
{
	glPushMatrix();
	glTranslatef(center.x(), center.y(), center.z());
	GLUquadricObj *quadObj = gluNewQuadric();

	gluQuadricDrawStyle(quadObj, GLU_FILL);
	gluQuadricNormals(quadObj, GLU_SMOOTH);

	gluSphere(quadObj, radius, 16, 16);

	gluDeleteQuadric(quadObj);
	glPopMatrix();
}

void SimpleDraw::DrawSpheres( StdVector<Vec3d> & centers, float radius /*= 1.0*/ )
{
	GLUquadricObj *quadObj = gluNewQuadric();

	gluQuadricDrawStyle(quadObj, GLU_FILL);
	gluQuadricNormals(quadObj, GLU_SMOOTH);

	for(int i = 0; i < (int)centers.size(); i++)
	{
		glPushMatrix();
		glTranslatef(centers[i].x(), centers[i].y(), centers[i].z());
		gluSphere(quadObj, radius, 16, 16);
		glPopMatrix();
	}

	gluDeleteQuadric(quadObj);
}

void SimpleDraw::DrawCylinder( const Vec3d & center, const Vec3d & direction /*= Vec3d (0,0,1)*/,
							  float height, float radius /*= 1.0f*/, float radius2 /*= -1*/ )
{
	glPushMatrix();
	glTranslatef(center.x(), center.y(), center.z());
	glMultMatrixd(qglviewer::Quaternion(qglviewer::Vec(0,0,1), qglviewer::Vec(direction)).matrix());

	GLUquadricObj *quadObj = gluNewQuadric();
	gluQuadricDrawStyle(quadObj, GLU_FILL);
	gluQuadricNormals(quadObj, GLU_SMOOTH);

	if(radius2 < 0)	radius2 = radius;

	gluCylinder(quadObj, radius, radius2, height, 16, 16);

	gluDeleteQuadric(quadObj);
	glPopMatrix();
}

void SimpleDraw::DrawArrow( Vec3d  from, Vec3d  to, bool isForward /*= true*/ , bool isFilledBase)
{
	if(!isForward){
		Vec3d  temp = from;
		from = to;
		to = temp;
	}

	if(isFilledBase) 
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	else
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	float length = (from-to).norm();
	float radius = length * 0.02f;
	if (radius < 0.0)
		radius = 0.05 * length;

	// draw cube base
	DrawCube(from, radius * 3);

	glPushMatrix();
	glTranslatef(from[0],from[1],from[2]);
	glMultMatrixd(qglviewer::Quaternion(qglviewer::Vec(0,0,1), qglviewer::Vec(to-from)).matrix());

	static GLUquadric* quadric = gluNewQuadric();

	const float head = 2.5*(radius / length) + 0.1;
	const float coneRadiusCoef = 4.0 - 5.0 * head;

	gluCylinder(quadric, radius, radius, length * (1.0 - head/coneRadiusCoef), 16, 1);
	glTranslatef(0.0, 0.0, length * (1.0 - head));
	gluCylinder(quadric, coneRadiusCoef * radius, 0.0, head * length, 16, 1);

	glPopMatrix();
}

void SimpleDraw::DrawArrowDirected( const Vec3d & pos, const Vec3d & normal, float height /*= 1.0f*/, 
								   bool isForward /*= true*/ , bool isFilledBase)
{
	DrawArrow(pos, pos + (normal*height), isForward, isFilledBase);
}

void SimpleDraw::PointArrowAt( Vec3d  point, float radius /*= 1.0f*/ )
{
	DrawArrowDirected(point, point + (double(radius) * Vec3d (Max(0.2, point.x()),point.y(),point.z()).normalized()), radius, false);
}

void SimpleDraw::DrawArrowDoubleDirected( const Vec3d & pos, const Vec3d & normal, float height /*= 1.0f*/, 
										 bool isForward /*= true*/, bool isFilledBase /*= true*/ )
{
	DrawArrowDirected(pos, normal, height, isForward, isFilledBase);

	glColor3f(1,0,0);
	DrawArrowDirected(pos, -normal, height, isForward, isFilledBase);
}

/*void SimpleDraw::IdentifyLines(StdVector<Line> & lines, float lineWidth, float r, float g, float b)
{
	glDisable(GL_LIGHTING);
	glLineWidth(lineWidth);

	glColor3f(r, g, b);

	glBegin(GL_LINES);
	for(int i = 0; i < (int)lines.size(); i++)
	{
		glVertex3dv(lines[i].a);
		glVertex3dv(lines[i].b);
	}
	glEnd();

	glEnable(GL_LIGHTING);
}*/
void SimpleDraw::IdentifyLine( const Vec3d & p1, const Vec3d & p2, Vec4d c, bool showVec3ds /*= true*/, float lineWidth /*= 3.0f*/ )
{
	glDisable(GL_LIGHTING);

	// Set color
	glColor4dv(c);

	glLineWidth(lineWidth);
	glBegin(GL_LINES);
	glVertex3dv(p1);
	glVertex3dv(p2);
	glEnd();

	if(showVec3ds)
	{
		// Draw colored end points
		glPointSize(lineWidth * 5);
		glBegin(GL_POINTS);
		glVertex3dv(p1);
		glVertex3dv(p2);
		glEnd();

		// White border end points
		glPointSize((lineWidth * 5) + 2);
		glColor3f(1, 1, 1);

		glBegin(GL_POINTS);
		glVertex3dv(p1);
		glVertex3dv(p2);
		glEnd();
	}

	glEnable(GL_LIGHTING);
}

void SimpleDraw::IdentifyLines( const StdVector<Vec3d> & p1, const StdVector<Vec3d> & p2, Vec4d c /*= Vec4d(1,0,0,1)*/, bool showVec3ds /*= true*/, float lineWidth /*= 3.0f */ )
{
	glDisable(GL_LIGHTING);

	// Set color
	glColor4dv(c);

	glLineWidth(lineWidth);
	glBegin(GL_LINES);
	for(int i = 0; i < (int)p1.size(); i++){
		glVertex3dv(p1[i]);
		glVertex3dv(p2[i]);
	}
	glEnd();

	if(showVec3ds)
	{
		// Draw colored end points
		glPointSize(lineWidth * 5);
		glBegin(GL_POINTS);
		for(int i = 0; i < (int)p1.size(); i++){
			glVertex3dv(p1[i]);
			glVertex3dv(p2[i]);
		}
		glEnd();

		// White border end points
		glPointSize((lineWidth * 5) + 2);
		glColor3f(1, 1, 1);

		glBegin(GL_POINTS);
		for(int i = 0; i < (int)p1.size(); i++){
			glVertex3dv(p1[i]);
			glVertex3dv(p2[i]);
		}
		glEnd();
	}

	glEnable(GL_LIGHTING);
}


void SimpleDraw::IdentifyPoint( const Vec3d & p, float r /*= 1.0*/, float g /*= 0.2f*/, float b /*= 0.2f*/, float pointSize /*= 10.0*/ )
{
	glDisable(GL_LIGHTING);

	// Colored dot
	glColor3f(r, g, b);
	glPointSize(pointSize);
	glBegin(GL_POINTS);
	glVertex3f(p.x(), p.y(), p.z());
	glEnd();

	// White Border
	glPointSize(pointSize + 2);
	glColor3f(1, 1, 1);

	glBegin(GL_POINTS);
	glVertex3f(p.x(), p.y(), p.z());
	glEnd();

	glEnable(GL_LIGHTING);
}

void SimpleDraw::IdentifyPoint2( Vec3d  p )
{
	// Green
	IdentifyPoint(p, 0.2f, 1.0f, 0.2f, 12.0f);
}

void SimpleDraw::IdentifyPoints(const StdVector<Vec3d > & points, Vec4d c, float pointSize)
{
	glDisable(GL_LIGHTING);

	glEnable(GL_BLEND); 
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// Colored dot
	glColor4dv(c);
	glPointSize(pointSize);
	glBegin(GL_POINTS);
	for(unsigned int i = 0; i < points.size(); i++)
		glVertex3dv(points[i]);
	glEnd();

	// White Border
	glPointSize(pointSize + 2);
	glColor4d(1, 1, 1, c[3]);

	glBegin(GL_POINTS);
	for(unsigned int i = 0; i < points.size(); i++)
		glVertex3dv(points[i]);
	glEnd();

	glEnable(GL_LIGHTING);
}

void SimpleDraw::IdentifyCurve( StdVector<Vec3d> & points, float r, float g, float b, float a, float lineWidth )
{
	glDisable(GL_LIGHTING);
	glColor4f(r, g, b, a);

	glLineWidth(lineWidth);

	glBegin(GL_LINE_STRIP);
	foreach(Vec3d p, points)
		glVertex3dv(p);
	glEnd();
}

void SimpleDraw::IdentifyConnectedPoints( StdVector<Vec3d> & points, float r /*= 0.4f*/, float g /*= 1.0*/, float b /*= 0.2f*/ )
{
	glDisable(GL_LIGHTING);
	glColor4f(r, g, b, 1);

	glLineWidth(2.0 + (r * 2) + (g * 4) + (b * 6));

	glBegin(GL_LINE_STRIP);
		foreach(Vec3d p, points)
			glVertex3dv(p);
	glEnd();

}

void SimpleDraw::IdentifyConnectedPoints2( StdVector<Vec3d > & points, float r /*= 0.4f*/, float g /*= 1.0*/, float b /*= 0.2f*/ )
{
	glDisable(GL_LIGHTING);

	int N = points.size();

	glLineWidth(3.0);
	glBegin(GL_LINE_STRIP);
	for(int i = 0; i < N; i++)
	{
		float t = Min((float(i) / N + 0.25f), 1.0);
		glColor3f(r * t, g * t, b * t);
		glVertex3dv(points[i]);
	}
	glEnd();

	// Colored dot
	glColor3f(r, g, b);
	glPointSize(13.0);
	glBegin(GL_POINTS);
	for(int i = 0; i < N; i++)
	{
		float t = float(i) / N;
		glColor3f(r * t, g * t, b * t);
		glVertex3dv(points[i]);
	}
	glEnd();

	// White Border
	glPointSize(15.0);
	glColor3f(1, 1, 1);

	glBegin(GL_POINTS);
	for(int i = 0; i < N; i++)
		glVertex3dv(points[i]);
	glEnd();

	glEnable(GL_LIGHTING);
}

void SimpleDraw::IdentifyLineRed( const Vec3d & p1, const Vec3d & p2, bool showVec3ds /*= true*/ )
{
	// Red line
	IdentifyLine(p1, p2, Vec4d(1.0, 0.2, 0.2, 1), showVec3ds);
}

void SimpleDraw::IdentifyArrow( const Vec3d  & start, const Vec3d  & end, float lineWidth /*= 2.0*/, float r /*= 1.0*/, float g /*= 0.2f*/, float b /*= 0.2f*/ )
{
	glDisable(GL_LIGHTING);

	// Transparency
	glEnable(GL_BLEND); 
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glLineWidth(lineWidth);

	glBegin(GL_LINES);

	glColor4f(r/2, g/2, b/2, 0.2f);
	glVertex3dv(start);

	glColor4f(r, g, b, 1.0);
	glVertex3dv(end);

	glEnd();

	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glLineWidth(1.0);
}

void SimpleDraw::IdentifyArrows( StdVector<Vec3d > & starts, StdVector<Vec3d > & ends, float lineWidth /*= 2.0*/, float r /*= 1.0*/, float g /*= 0.2f*/, float b /*= 0.2f*/ )
{
	glDisable(GL_LIGHTING);

	// Transparency
	glEnable(GL_BLEND); 
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glLineWidth(lineWidth);

	glBegin(GL_LINES);
	for(unsigned int i = 0; i < starts.size(); i++)
	{
		glColor4f(r, g, b, 0.0);
		glVertex3dv(starts[i]);

		glColor4f(r, g, b, 1.0);
		glVertex3dv(ends[i]);
	}
	glEnd();

	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
}

void SimpleDraw::DrawBarChart(const StdVector<double> & data, int x, int y, double height, double width, int barWidth)
{
	glDisable(GL_LIGHTING);
	glLineWidth(barWidth);

	glPushMatrix();
	glTranslated(x,y,0);

	double scalingY = height / (1 + *max_element(data.begin(), data.end()));
	double scalingX = width * (1.0 / (int)data.size());

	glBegin(GL_LINES);
	for(int i = 0; i < (int)data.size(); i++)
	{
		glVertex2i(i * barWidth * scalingX, 0);
		glVertex2i(i * barWidth * scalingX, -data[i] * scalingY);
	}
	glEnd();

	glPopMatrix();

	glEnable(GL_LIGHTING);
}

void SimpleDraw::drawCornerAxis(const double * cameraOrientation)
{
	int viewport[4];
	int scissor[4];

	// The viewport and the scissor are changed to fit the lower left
	// corner. Original values are saved.
	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetIntegerv(GL_SCISSOR_BOX, scissor);

	// Axis viewport size, in pixels
	const int size = 50;
	glViewport(0,0,size,size);
	glScissor(0,0,size,size);

	// The Z-buffer is cleared to make the axis appear over the
	// original image.
	glClear(GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(-1, 1, -1, 1, -1, 1);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glMultMatrixd(cameraOrientation);

	// Tune for best line rendering
	glDisable(GL_LIGHTING);
	glLineWidth(3.0);

	glBegin(GL_LINES);
	glColor3f(1.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(1.0, 0.0, 0.0);

	glColor3f(0.0, 1.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 1.0, 0.0);

	glColor3f(0.0, 0.0, 1.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, 1.0);
	glEnd();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glEnable(GL_LIGHTING);

	// The viewport and the scissor are restored.
	glScissor(scissor[0],scissor[1],scissor[2],scissor[3]);
	glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);
}

std::vector<Vec4d> SimpleDraw::RandomColors( int count )
{
	std::vector<Vec4d> colors;
	for (int i = 0; i < count; i++){
		float r = ((rand() % 225) + 30) / 255.0f;
		float g = ((rand() % 230) + 25) / 255.0f;
		float b = ((rand() % 235) + 20) / 255.0f;
		colors.push_back(Vec4d(r, g, b, 1.0));
	}
	return colors;
}

void SimpleDraw::DrawGraph2D( const StdVector< StdVector <double> > & data, double maximum)
{
	int w = data.front().size();
	int h = data.size();

	// Compute length of each grid quad
	double q = 1.0 / Max(w,h);

	// Find actual data max
	double actualMax = DBL_MIN;
	for(int i = 0; i < (int)data.size(); i++)
		actualMax = Max(actualMax, MaxElement(data[i]));

	uchar rgb[3];	

	glDisable(GL_LIGHTING);

	glBegin(GL_QUADS);

	for(int y = 0; y < h - 1; y++)
	{
		for(int x = 0; x < w - 1; x++)
		{
			double c1 = Max(0, data[y][x]);
			double c2 = Max(0, data[y+1][x]);
			double c3 = Max(0, data[y][x+1]);
			double c4 = Max(0, data[y+1][x+1]);

			ColorMap::jetColorMap(rgb, c1, 0, actualMax);
			glColor3ubv(rgb);
			glVertex3dv(Vec3d(y * q, x * q, c1 / maximum));

			ColorMap::jetColorMap(rgb, c3, 0, actualMax);
			glColor3ubv(rgb);
			glVertex3dv(Vec3d(y * q, (x+1) * q, c3 / maximum));

			ColorMap::jetColorMap(rgb, c4, 0, actualMax);
			glColor3ubv(rgb);
			glVertex3dv(Vec3d((y+1) * q, (x+1) * q, c4 / maximum));

			ColorMap::jetColorMap(rgb, c2, 0, actualMax);
			glColor3ubv(rgb);
			glVertex3dv(Vec3d((y+1) * q, x * q, c2 / maximum));
		}
	}

	glEnd();

	glEnable(GL_LIGHTING);
}

void SimpleDraw::DrawCircle( const Vec3d& center, double radius, const Vec4d& c, const Vec3d& n, float lineWidth )
{
	Vec3d startV(0,0,0);

	// Find orthogonal start vector
	if ((abs(n.y()) >= 0.9f * abs(n.x())) && 
		abs(n.z()) >= 0.9f * abs(n.x())) startV = Vec3d(0.0f, -n.z(), n.y());
	else if ( abs(n.x()) >= 0.9f * abs(n.y()) && 
		abs(n.z()) >= 0.9f * abs(n.y()) ) startV = Vec3d(-n.z(), 0.0f, n.x());
	else startV = Vec3d(-n.y(), n.x(), 0.0f);

	int segCount = 20;
	double theta = 2.0 * M_PI / segCount;

	glDisable(GL_LIGHTING);
	glLineWidth(lineWidth);
	glColor4dv(c);

	glBegin(GL_LINE_LOOP);
	for(int i = 0; i < segCount; i++){
		glVertex3dv(center + startV * radius );
		ROTATE_VEC(startV, theta, n);
	}
	glEnd();

	glEnable(GL_LIGHTING);
}

void SimpleDraw::DrawCircles( const StdVector<Vec3d>& centers, const StdVector<double>& radius, const Vec4d& c, const Vec3d& normal, float lineWidth )
{
	for(int i = 0; i < (int) centers.size(); i++)
		DrawCircle(centers[i], radius[i], c, normal, lineWidth);
}
