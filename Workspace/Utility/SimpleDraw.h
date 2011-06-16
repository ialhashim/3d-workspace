#ifndef SIMPLE_DRAW_H
#define SIMPLE_DRAW_H

#include "Line.h"

#ifndef Vector
#include <vector>
#define Vector std::vector
#endif

class SimpleDraw
{
public:

	// POINTS
	static void IdentifyPoint(const Vec& p, float r = 1.0, float g = 0.2f, float b = 0.2f, float pointSize = 10.0);
	static void IdentifyPoint2(Vec p);
	static void IdentifyPoints(Vector<Vec> & points, float r = 1.0, float g = 0.2f, float b = 0.2f, float pointSize = 10.0);
	static void IdentifyConnectedPoints(Vector<Vec> & points, float r = 0.4f, float g = 1.0, float b = 0.2f);

	// LINES
	static void IdentifyLine(const Vec& p1, const Vec& p2, bool showPoints = true);
	static void IdentifyLineRed(const Vec& p1, const Vec& p2, bool showPoints = true);
	static void IdentifyLine(const Vec& p1, const Vec& p2, float r, float g, float b, bool showPoints = true, float lineWidth = 3.0f);
	static void IdentifyLines(Vector<Line> & lines, float lineWidth = 1.0, float r = 1.0, float g = 0.6f, float b = 0);

	// Primitives
	static void DrawCube(const Vec& center, float length = 1.0);
	static void DrawSphere(const Vec& center, float radius = 1.0);
	static void DrawCylinder(const Vec& center, const Vec& direction = Vec(0,0,1), float height = 1.0, float radius = 1.0, float radius2 = -1);

	static void DrawBox(const Vec& center, float width, float length, float height, float r = 0, float g = 1.0, float b = 0);
	static void DrawSolidBox(const Vec& center, float width, float length, float height, float r = 0, float g = 1.0, float b = 0, float a = 1.0);
	static void DrawTriangle(const Vec& v1, const Vec& v2, const Vec& v3,float r = 0.2f, float g = 1.0, float b = 0.1f, float a = 1.0, bool isOpaque = true);
	static void DrawTriangles(const Vector< Vector<Vec> > & tris, float r = 0.2f, float g = 1.0, float b = 0.1f, float a = 1.0, bool isOpaque = true, bool isDrawPoints = true);

	static void DrawSquare(const Vec& v1, const Vec& v2, const Vec& v3, const Vec& v4, bool isOpaque = true, float r = 0.1f, float g = 0.2f, float b = 1.0, float a = 1.0);
	static void DrawSquares(const Vector<Vector<Vec> > & squares, bool isOpaque = true, float r = 0.1f, float g = 0.2f, float b = 1.0, float a = 1.0);
	static void DrawLineTick(const Vector<Vec>& start, const  Vector<Vec>& direction, float len = 0.25f, bool border = false, float r = 0.65f, float g = 0.6f, float b = 0.8f, float a = 0.9f);

	// ARROWS
	static void IdentifyArrow(const Vec & start, const Vec & end, float lineWidth = 2.0, float r = 1.0, float g = 0.2f, float b = 0.2f);
	static void IdentifyArrows(Vector<Vec> & starts, Vector<Vec> & ends, float lineWidth = 2.0, float r = 1.0, float g = 0.2f, float b = 0.2f);

	static void DrawArrow(Vec from, Vec to, bool isForward = true, bool isFilledBase = true);
	static void DrawArrowDirected(const Vec& pos, const Vec& normal, float height = 1.0, bool isForward = true, bool isFilledBase = true);
	static void DrawArrowDoubleDirected(const Vec& pos, const Vec& normal, float height = 1.0, bool isForward = true, bool isFilledBase = true);

	static void PointArrowAt(Vec point, float radius = 1.0);
	static void DrawBarChart(const Vector<double> & data, int x, int y, double height = 100.0, double width = 300, int barWidth = 1.0);
	
	// Misc.
	static void drawCornerAxis(const double * cameraOrientation);
};

#endif // SIMPLE_DRAW_H
