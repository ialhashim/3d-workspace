#include "Spline.h"

Spline::Spline(): vp(), delta_t(0)
{
}

Spline::Spline(const Spline& s)
{
	for (int i = 0; i < (int)s.vp.size(); i++)
		vp.push_back(s.vp[i]);
	delta_t = s.delta_t;
}

Spline& Spline::operator= (const Spline& s)
{
	for (int i = 0; i < (int)s.vp.size(); i++)
		this->vp.push_back(s.vp[i]);
	this->delta_t = s.delta_t;

	return *this;
}

// Solve the Catmull-Rom parametric equation for a given time(t) and vector quadruple (p1,p2,p3,p4)
Vec Spline::Eq(double t, const Vec& p1, const Vec& p2, const Vec& p3, const Vec& p4)
{
	double t2 = t * t;
	double t3 = t2 * t;

	double b1 = .5 * (  -t3 + 2*t2 - t);
	double b2 = .5 * ( 3*t3 - 5*t2 + 2);
	double b3 = .5 * (-3*t3 + 4*t2 + t);
	double b4 = .5 * (   t3 -   t2    );

	return (p1*b1 + p2*b2 + p3*b3 + p4*b4); 
}

void Spline::AddSplinePoint(const Vec& v)
{
	vp.push_back(v);
	delta_t = (double)1 / (double)vp.size();
}

Vec Spline::GetInterpolatedSplinePoint(double t)
{
	// Find out in which interval we are on the spline
	int p = (int)(t / delta_t);
	// Compute local control point indices
#define BOUNDS(pp) { if (pp < 0) pp = 0; else if (pp >= (int)vp.size()-1) pp = vp.size() - 1; }
	int p0 = p - 1;     BOUNDS(p0);
	int p1 = p;         BOUNDS(p1);
	int p2 = p + 1;     BOUNDS(p2);
	int p3 = p + 2;     BOUNDS(p3);
	// Relative (local) time 
	double lt = (t - delta_t*(double)p) / delta_t;
	// Interpolate
	return Spline::Eq(lt, vp[p0], vp[p1], vp[p2], vp[p3]);
}

int Spline::GetNumPoints()
{
	return vp.size();
}

Vec& Spline::GetNthPoint(int n)
{
	return vp[n];
}

double* Spline::GetLastPoint()
{
	return &vp[vp.size() - 1].x;
}

Vector<Vec> Spline::GetPoints()
{
	return vp;
}

void Spline::Translate(const Vec & to)
{
	Vec delta = to - vp[0];

	for(int i = 0; i < (int)vp.size(); i++)
		vp[i] += delta;
}

void Spline::Align(const Vec & direction)
{
	qglviewer::Quaternion q((vp[1] - vp[0]).unit(), direction.unit());

	for(int i = 0; i < (int)vp.size(); i++)
		vp[i] = q.rotate(vp[i]);
}

void Spline::Crop(int start, int end)
{
	Vector<Vec> temp_vp;

	temp_vp.push_back(vp[start]);

	for(int i = start; i < (int)vp.size() - end; i++)
		temp_vp.push_back(vp[i]);

	temp_vp.push_back(vp[(vp.size() - end) - 1]);

	vp = temp_vp;

	delta_t = (double)1 / (double)vp.size();
}

void Spline::Subdivide(int numIterations)
{
	for(int s = 0; s < numIterations; s++)
	{
		Vector<Vec> temp_vp;

		// First
		temp_vp.push_back(vp[0]);

		for(int i = 0; i < (int)vp.size() - 1; i++)
		{
			temp_vp.push_back( Point3D::MidVec(vp[i], vp[i+1]) );
			temp_vp.push_back( vp[i+1] );
		}

		// Last
		temp_vp.push_back(vp[vp.size() - 1]);

		vp = temp_vp;

		Smooth(4);
	}

	delta_t = (double)1 / (double)vp.size();
}

void Spline::Simplify(int numIterations)
{
	for(int s = 0; s < numIterations; s++)
	{
		Vector<Vec> temp_vp;

		if(vp.size() < 4)
			return;

		for(int i = 0; i < (int)vp.size(); i += 2)
			temp_vp.push_back(vp[i]);

		vp = temp_vp;
	}

	delta_t = (double)1 / (double)vp.size();
}

void Spline::Smooth(int numIterations)
{
	for(int s = 0; s < numIterations; s++)
	{
		Vector<Vec> temp_vp;

		for(int i = 0; i < (int)vp.size(); i++)
		{
			if(i == 0 || i == (int)vp.size() - 1)
				temp_vp.push_back(vp[i]);
			else
				temp_vp.push_back( (vp[i - 1] + vp[i + 1]) / 2.0f);
		}

		vp = temp_vp;
	}

	delta_t = (double)1 / (double)vp.size();
}

Vector<Vec> Spline::GetUniformPath(double stepLength, const Vector<Vec>& fromPoints)
{
	Vector<Vec> path;

	// Setup spline
	BezierSpline uniform_spline(fromPoints);

	// Make uniform
	BezierSpline temp = uniform_spline;
	uniform_spline.SetUniformDistanceSampled(temp, fromPoints.size());

	// Get path
	double dt = (double)stepLength / uniform_spline.GetTotalLength();

	for(double t = 0; t <= 1.0; t += dt)
		path.push_back(uniform_spline.GetValue( t ));

	return path;
}

Vector<Vec> Spline::GetUniformPath( int numSteps, const Vector<Vec>& fromPoints )
{
	Vector<Vec> path;

	// Setup spline
	BezierSpline uniform_spline(fromPoints);

	// Make uniform
	BezierSpline temp = uniform_spline;
	uniform_spline.SetUniformDistanceSampled(temp, fromPoints.size());

	double delta = 1.0f / numSteps;

	for(double i = 0; i <= 1.0f; i += delta)
	{
		path.push_back(uniform_spline.GetValue( i ) );
	}

	return path;
}

Vector<Vec> Spline::GetUniformPathSteps(int numSteps)
{
	Vector<Vec> path;

	// Setup spline
	BezierSpline uniform_spline(vp, true);

	// Make uniform
	BezierSpline temp = uniform_spline;
	uniform_spline.SetUniformDistanceSampled(temp, vp.size());

	double delta = 1.0f / numSteps;

	for(double i = 0; i <= 1.0f; i += delta)
	{
		path.push_back(uniform_spline.GetValue( i ) );
	}

	return path;
}

Vector<Vec> Spline::GetUniformPath(double stepLength)
{
	Vector<Vec> path;

	// Setup spline
	BezierSpline uniform_spline(vp);

	// Make uniform
	BezierSpline temp = uniform_spline;
	uniform_spline.SetUniformDistanceSampled(temp, vp.size());

	double dt = stepLength / uniform_spline.GetTotalLength();

	for(double t = 0; t <= 1; t += dt)
	{
		path.push_back(uniform_spline.GetValue( t ));
	}

	return path;
}

Vector<Vec> Spline::GetUniformTangent(double stepLength)
{
	Vector<Vec> tangent;

	// Setup spline
	BezierSpline uniform_spline(vp);

	// Make uniform
	BezierSpline temp = uniform_spline;
	uniform_spline.SetUniformDistanceSampled(temp, vp.size());

	// Get path
	double dt = stepLength / uniform_spline.GetTotalLength();

	for(double t = 0; t <= 1; t += dt)
	{
		tangent.push_back( uniform_spline.GetDerivativeGlobal(t).unit() );
	}

	return tangent;
}

Vector<Vec> Spline::GetUniformTangent( int numSteps, const Vector<Vec>& fromPoints )
{
	Vector<Vec> tangent;

	// Setup spline
	BezierSpline uniform_spline(fromPoints);

	// Make uniform
	BezierSpline temp = uniform_spline;
	uniform_spline.SetUniformDistanceSampled(temp, fromPoints.size());

	double delta = 1.0f / numSteps;

	for(double t = 0; t <= 1.0f; t += delta)
	{
		tangent.push_back( uniform_spline.GetDerivativeGlobal( t ).unit() );
	}

	return tangent;
}

void Spline::Clear()
{
	vp.clear();
	delta_t = 0;
}

BezierSpline Spline::ToUniformBezierSpline()
{
	// Setup spline
	BezierSpline uniform_spline(vp);

	// Make uniform
	BezierSpline temp = uniform_spline;
	uniform_spline.SetUniformDistanceSampled(temp, vp.size());

	return uniform_spline;
}
