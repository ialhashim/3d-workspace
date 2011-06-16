#pragma once

#include "Bezier.h"

// Based on Charles Bloom Galaxy3 http://www.cbloom.com/3d/galaxy3/index.html

enum EEndHandling{ eEnd_Clamp,	eEnd_Wrap, };
struct Knot{
	Vec	point;
	double	time;
};

struct PlaceOnSpline{
	int		segment;
	double	localTime;
};

struct Piece{
	Bezier		curve;
	double		accumTime;	// parameters : // time of the end of this piece
	double		accumLength;  // derived : // length + length of predecessors
};

enum ESplineType{eInvalid,eGeneric,eUniformTime,eUniformLength};

/*
	a BezierSpline is just a bunch of Bezier curves,
	with piecewise linear time control for how they fit together
*/

class BezierSpline
{
private:
	Vector<Piece>	m_pieces;
	ESplineType		m_type;

public:
	 BezierSpline();
	~BezierSpline();

	BezierSpline(const Vector<Vec> & pKnots, bool loop = false) 
	{ 
		if(loop)
			SetFromPoints(pKnots, eEnd_Wrap); 
		else
			SetFromPoints(pKnots);
	};

	bool IsValid() const;

	//----------------------------------------------
	// Setup :

	void SetFromPoints(const Vector<Vec> & pKnots, const EEndHandling eEnds = eEnd_Clamp);
	void SetFromKnots(const Vector<Knot> & pKnots, const EEndHandling eEnds = eEnd_Clamp);
	
	//! Sample another spline at equi-distant spots to make a new spline with
	//	equal arc-length spaced knots 
	//void SetUniformSampled(const BezierSpline & from,const int numKnots);

	void SetUniformTimeSampled(const BezierSpline & from,const int numKnots);
	void SetUniformDistanceSampled(const BezierSpline & from,const int numKnots);

	void SetSubdivisionSampled(const BezierSpline & from,const double maxErrSqr);
	void MakeSimplification(const BezierSpline & from,const double maxErrSqr);
	void MakeMinimumUniformTimeSampled(const BezierSpline & from,const double maxErrSqr);

	//----------------------------------------------
	// Queries :

	int GetNumSegments() const { return m_pieces.size(); }
	const Bezier & GetCurve(int seg) const { return m_pieces[seg].curve; }

	double GetTotalLength() const;

	//! the Length functions here treat arc-length as if it
	//	varies linearly on each segment (an approximation)
	void GetSegmentTimeRange(int seg,double *pStart,double *pEnd) const;
	void GetSegmentLengthRange(int seg,double *pStart,double *pEnd) const;

	void GetPlaceByTime(const double t,PlaceOnSpline * pInto) const;
	void GetPlaceByLength(const double length,PlaceOnSpline * pInto) const;
	void GetPlaceByLengthAccurate(const double length,PlaceOnSpline * pInto) const;

	double GetTimeFromPlace(const PlaceOnSpline & place) const;
	double GetLengthFromPlace(const PlaceOnSpline & place) const;

	double GetTimeClosestPoint(const Vec & point) const;
	Vec ClosestPoint( const Vec & point ) const;

	//----------------------------------------------
	//! To do serious queries, get the Curve and then talk to it
	//	here are some little helpers :
	const Vec GetValue(const double t) const
	{
		PlaceOnSpline place;
		GetPlaceByTime(t,&place);
		return GetCurve(place.segment).GetValue(place.localTime);	
	}
	const Vec GetValue(const PlaceOnSpline & place) const
	{
		return GetCurve(place.segment).GetValue(place.localTime);	
	}
	//! this is the derivative with respect to local time;
	//	if the knots are not equally spaced in time, this needs to be adjusted
	//	to be the correct derivative WRST global time
	const Vec GetDerivativeLocal(const double t) const
	{
		PlaceOnSpline place;
		GetPlaceByTime(t,&place);
		return GetCurve(place.segment).GetDerivative(place.localTime);	
	}
	const Vec GetDerivativeLocal(const PlaceOnSpline & place) const
	{
		return GetCurve(place.segment).GetDerivative(place.localTime);	
	}
	const Vec GetDerivativeGlobal(const double t) const
	{
		PlaceOnSpline place;
		GetPlaceByTime(t,&place);
		return GetDerivativeGlobal(place);
	}
	const Vec GetDerivativeGlobal(const PlaceOnSpline & place) const;

	static double DifferenceSqr(const BezierSpline & s1,const BezierSpline & s2,const int numSamples);

private:

	void SetSubdivisionSampled_RecursiveAdd( const Bezier & curve, const double t0, const double t1, const double maxErrSqr);

	void UpdateDerived();
	void ParameterizeByArcLength();
};

int IndexOf(const int i, const int size, const EEndHandling eEnds);
