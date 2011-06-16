#include "Bezier.h"

/*
	Controls are B0,B1,B2,B3. The Curve function is :

	s = (1-t)
	B(t) = B0 * s^3 + B1 * s^2 * t * 3 + B2 * s * t^2 * 3 + B3 * t^3
	B'(t) = 3 * [ (B1-B0) * s^2 + (B2-B1) * st + (B3 - B2) * t^2 ]
*/

#define B0	m_controls[0]
#define B1	m_controls[1]
#define B2	m_controls[2]
#define B3	m_controls[3]

void Bezier::SetFromControls(const Vec * pControls)
{
	memcpy(m_controls,pControls,4*sizeof(Vec));
}

void Bezier::SetFromControls(const Vec & b0,const Vec & b1,const Vec & b2,const Vec & b3)
{
	B0 = b0;
	B1 = b1;
	B2 = b2;
	B3 = b3;
}

void Bezier::SetFromEnds(const Vec & v0,const Vec & d0, const Vec & v1,const Vec & d1)
{
	B0 = v0;
	B3 = v1;
	//d0 = 3*(B1-B0);
	B1 = d0*(1.0/3.0) + B0;
	B2 = B3 - d1*(1.0/3.0);
}

const Vec Bezier::GetDerivative0() const
{
	return 3.0*(B1 - B0);
}

const Vec Bezier::GetDerivative1() const
{
	return 3.0*(B3 - B2);
}

const Vec Bezier::Get2ndDerivative0() const
{
	return 6.0*(B2 - B0);
}

const Vec Bezier::Get2ndDerivative1() const
{
	return 6.0*(B3 - B1);
}

const Vec Bezier::Get3rdDerivative() const
{
	return 6.0*( B0 - B1 - B2 + B3 );
}

const Vec Bezier::GetValue(const double t) const
{
	const double s = (1.0-t);
	return B0* (s*s*s) + B1*(s*s*t*3.0) + B2*(t*t*s*3.0) + B3*(t*t*t);
}

const Vec Bezier::GetDerivative(const double t) const
{
	const double s = (1.0-t);
	return (B1-B0) * (3.0*s*s) + (B2-B1) * (6.0*s*t) + (B3 - B2) * (3.0*t*t);
}

const Vec Bezier::Get2ndDerivative(const double t) const
{
	return MakeLerp( Get2ndDerivative0(),Get2ndDerivative1(), t );
}

double Bezier::ComputePartialLength(const double start,const double end,const double dt /*= 0.01f*/) const
{
	// Length requires an integration!!
	assert( fiszerotoone(start) );
	assert( fiszerotoone(end) );

	Vec pos[2];
	int iPos=1;
	pos[0] = GetValue(start);

	double length = 0.0;

	double t;
	for(t = start + dt; t <= end; t += dt)
	{
		// sample at t into pos[iPos];
		pos[iPos] = GetValue(t);
		// add on current distance
		const double d = Distance(pos[0], pos[1]);
		length += d;
		// swap iPos
		iPos ^= 1;
	}

	// add on the last bit :
	if ( t < end )
	{
		pos[iPos] = GetValue(end);
		const double d = Distance(pos[0], pos[1]);
		length += d;		
	}

	return length;
}

double Bezier::ComputeParameterForLength(const double lengthParam,const double dt /*= 0.01f*/) const
{
	if ( lengthParam <= 0.0 )
		return 0.0;

	Vec pos[2];
	int iPos=1;
	pos[0] = GetValue(0.0);

	double length = 0.0;

	for(double t = 0.0; t <= 1.0; t += dt)
	{
		// sample at t into pos[iPos];
		pos[iPos] = GetValue(t+dt);
		// add on current distance
		const double d = Distance(pos[0],pos[1]);
		const double nextLength = length + d;

		if ( lengthParam <= nextLength )
		{
			// I'm somewhere in between (t) and (t+dt)
			double lerper = fmakelerperclamp(length,nextLength,lengthParam);
			return flerp(t,t+dt,lerper);
		}

		length = nextLength;
		// swap iPos
		iPos ^= 1;
	}

	return 1.0;
}

void Bezier::Subdivide(Bezier *pBez1,Bezier *pBez2) const
{
	Vec b01,b12,b23,b012,b123,b0123;

	SetAverage(b01,B0,B1);
	SetAverage(b12,B1,B2);
	SetAverage(b23,B2,B3);
	SetAverage(b012,b01,b12);
	SetAverage(b123,b12,b23);
	SetAverage(b0123,b012,b123);
	
	pBez1->SetFromControls(B0,b01,b012,b0123);
	pBez2->SetFromControls(b0123,b123,b23,B3);
}

void Bezier::Subdivide(Bezier *pBez1,Bezier *pBez2,const double t) const
{
	Vec b01,b12,b23,b012,b123,b0123;

	SetLerp(b01,B0,B1,t);
	SetLerp(b12,B1,B2,t);
	SetLerp(b23,B2,B3,t);
	SetLerp(b012,b01,b12,t);
	SetLerp(b123,b12,b23,t);
	SetLerp(b0123,b012,b123,t);
	
	pBez1->SetFromControls(B0,b01,b012,b0123);
	pBez2->SetFromControls(b0123,b123,b23,B3);
}

//! tells how far the Bezier differs from a line - in units squared
double Bezier::ComputeLinearErrorBoundSqr() const
{
	/*
		Error is bounded by
		[ B1 - lerp(B0,B3,(1/3)) ] or [ B2 - ..(2/3) ]

		Actually, there's a tighter bound, which is
		[ B1 - projection of B1 onto (B0,B3) ]
		but this is much cheaper
	*/

	const Vec d1 = B0 * (2.0/3.0) + B3 * (1.0/3.0) - B1;
	const Vec d2 = B0 * (1.0/3.0) + B3 * (2.0/3.0) - B2;

	const double errSqr = gMAX( d1.squaredNorm(), d2.squaredNorm() );

	return errSqr;
}

double Bezier::MinDistanceToControlPoints(const Vec & point) const
{
	double dist1 = (B0 - point).norm();
	double dist2 = (B1 - point).norm();
	double dist3 = (B2 - point).norm();
	double dist4 = (B3 - point).norm();

	return Min(dist1, Min(dist2, Min(dist3, dist4)));
}
