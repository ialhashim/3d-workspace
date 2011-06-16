#pragma once

#include "BezierUtil.h"

#define EPSIL 0.001

// Based on Charles Bloom Galaxy3 http://www.cbloom.com/3d/galaxy3/index.html

// Cubic Bezier curve parameterized in 0->1
class Bezier
{
public:
	 Bezier() { }
	~Bezier() { }

	void SetFromControls(const Vec * pControls);
	void SetFromControls(const Vec & b0,const Vec & b1,const Vec & b2,const Vec & b3);
	void SetFromEnds(const Vec & v0,const Vec & d0,const Vec & v1,const Vec & d1);

	const Vec * GetControls() const { return m_controls; }

	const Vec GetValue(const double t) const;
	const Vec GetDerivative(const double t) const;
	const Vec Get2ndDerivative(const double t) const;
	const Vec Get3rdDerivative() const;

	const Vec GetValue0() const { return m_controls[0]; }
	const Vec GetValue1() const { return m_controls[3]; }
	const Vec GetDerivative0() const;
	const Vec GetDerivative1() const;
	const Vec Get2ndDerivative0() const;
	const Vec Get2ndDerivative1() const;

	// length is computed with an iterative integration !
	double ComputePartialLength(const double start,const double end,const double dt = EPSIL) const;
	double ComputeLength(const double dt = EPSIL) const { return ComputePartialLength(0.0,1.0,dt); }
	double ComputeParameterForLength(const double length,const double dt = EPSIL) const;

	// subdivide makes two beziers for the [0,0.5] and [.5,1] ranges
	void Subdivide(Bezier *pBez1,Bezier *pBez2) const;
	
	// breaks [0,t][t,1]
	void Subdivide(Bezier *pBez1,Bezier *pBez2,const double t) const;

	//! tells how far the Bezier differs from a line - in units squared
	double ComputeLinearErrorBoundSqr() const;

	double MinDistanceToControlPoints(const Vec & point) const;

private:
	Vec	m_controls[4];
	// 3*4*4 = 48 bytes
};
