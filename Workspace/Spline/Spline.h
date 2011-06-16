#pragma once

#include "BezierSpline.h"

// Based on http://www.codeproject.com/KB/recipes/Overhauser.aspx
class Spline
{
private:
    Vector<Vec> vp;
    double delta_t;

public:

    // Constructors and destructor
    Spline();
    Spline(const Spline&);
	Spline& operator= (const Spline&); 

    // Operations
    void AddSplinePoint(const Vec& v);
	Vec GetInterpolatedSplinePoint(double t);   // t = 0...1; 0=vp[0] ... 1=vp[max]
	Vec& GetNthPoint(int n);
	double* GetLastPoint();
	Vector<Vec> GetPoints();

	int GetNumPoints();
	
    // Static method for computing the Catmull-Rom parametric equation
    // given a time (t) and a vector quadruple (p1,p2,p3,p4).
    static Vec Eq(double t, const Vec& p1, const Vec& p2, const Vec& p3, const Vec& p4);

	// Transformation
	void Translate(const Vec & to);
	void Align(const Vec & direction);

	// Modifications
	void Crop(int start, int end);
	void Simplify(int numIterations);
	void Subdivide(int numIterations);
	void Smooth(int numIterations);
	void Clear();

	// Generate paths & tangents
	Vector<Vec> GetUniformPathSteps(int numSteps);

	Vector<Vec> GetUniformPath(double stepLength);
	Vector<Vec> GetUniformTangent(double stepLength);

	BezierSpline ToUniformBezierSpline();

	// Static
	static Vector<Vec> GetUniformPath(double stepLength, const Vector<Vec>& fromPoints);
	static Vector<Vec> GetUniformPath(int numSteps, const Vector<Vec>& fromPoints);
	static Vector<Vec> GetUniformTangent( int numSteps, const Vector<Vec>& fromPoints );
};
