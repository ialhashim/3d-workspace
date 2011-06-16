#include "BezierSpline.h"

BezierSpline::BezierSpline() : m_type(eInvalid)
{
}

BezierSpline::~BezierSpline()
{
}

void BezierSpline::SetFromPoints(const Vector<Vec> & pKnots, const EEndHandling eEnds)
{
	int numKnots = pKnots.size();

	int numCurves = numKnots-1;
	assert( numCurves > 0 );

	// make curves from knots :
	m_pieces.clear();
	m_pieces.resize(numCurves);
	
	for(int i = 0; i < numCurves; i++)
	{
		// P0 = pKnots[i];
		// P1 = pKnots[i+1];
		// D0 = (pKnots[i+1]-pKnots[i-1])/6
		// D1 = (pKnots[i+2]-pKnots[i])/6

		Vec D0 = ( pKnots[i+1] - pKnots[IndexOf(i-1, numKnots, eEnds)] ) / 2.f;
		Vec D1 = ( pKnots[IndexOf(i+2, numKnots, eEnds)] - pKnots[i] ) / 2.f;

		m_pieces[i].curve.SetFromEnds(pKnots[i], D0, pKnots[i+1], D1);
		m_pieces[i].accumTime = (i + 1) / double(numCurves);
	}

	m_type = eUniformTime;

	UpdateDerived();
	// do NOT ParameterizeByLength, need even time params to match derivatives
}

//! SetFromKnots is nearly SetFromPoints, but the time parameter is provided
void BezierSpline::SetFromKnots(const Vector<Knot> & pKnots, const EEndHandling eEnds)
{
	int numKnots = pKnots.size();

	int numCurves = numKnots-1;
	assert( numCurves > 0 );

	// make curves from knots :
	m_pieces.clear();
	m_pieces.resize(numCurves);
	
	assert( pKnots[0].time == 0.f );
	assert( pKnots[numKnots-1].time == 1.f );

	for(int i = 0; i < numCurves; i++)
	{
		// P0 = pKnots[i];
		// P1 = pKnots[i+1];
		// D0 = (pKnots[i+1]-pKnots[i-1])/6
		// D1 = (pKnots[i+2]-pKnots[i])/6

		Vec D0 = ( pKnots[i+1].point - pKnots[IndexOf(i-1,numKnots,eEnds)].point )/2.f;
		Vec D1 = ( pKnots[IndexOf(i+2,numKnots,eEnds)].point - pKnots[i].point )/2.f;

		double timeScale = (pKnots[i+1].time - pKnots[i].time) / double(numCurves);
		// if time was even, timeScale would be 1.0
		// if timeScale is > 1.0 it means we move slower on this knot than we would think,
		//	so then we have to scale up D to match the neighbor
		D0 *= timeScale;
		D1 *= timeScale;

		m_pieces[i].curve.SetFromEnds(pKnots[i].point,D0, pKnots[i+1].point,D1);
		m_pieces[i].accumTime = pKnots[i+1].time;
	}
	
	m_type = eGeneric;

	UpdateDerived();
}

//! Sample another spline at equi-distant spots to make a new spline with
//	equal arc-length spaced knots 
void BezierSpline::SetUniformTimeSampled(const BezierSpline & from, const int numKnots)
{
	int numCurves = numKnots-1;
	assert( numCurves > 0 );

	double timePerSeg = 1.f / double(numCurves);
	
	// make curves from knots :
	m_pieces.clear();
	m_pieces.resize(numCurves);
	
	Vec P[2];
	Vec D[2];
	int index=1;
	P[0] = from.GetValue(0.f);
	D[0] = from.GetDerivativeGlobal(0.f);

	for(int i = 0; i < numCurves; i++)
	{
		double t = (i+1) * timePerSeg;
		if ( i == (numCurves-1) )
			t = 1.f;
		PlaceOnSpline place;
		from.GetPlaceByTime(t,&place);

		P[index] = from.GetValue(place);
		D[index] = from.GetDerivativeGlobal(place);

		const Vec & P0 = P[index^1];
		const Vec & P1 = P[index];

		const Vec D0 = D[index^1] * timePerSeg;
		const Vec D1 = D[index] * timePerSeg;

		m_pieces[i].curve.SetFromEnds(P0,D0, P1,D1);
		m_pieces[i].accumTime = t;

		index ^=1;
	}

	m_type = eUniformTime;

	UpdateDerived();
}


//! Sample another spline at equi-distant spots to make a new spline with
//	equal arc-length spaced knots 
void BezierSpline::SetUniformDistanceSampled(const BezierSpline & from,const int numKnots)
{
	int numCurves = numKnots-1;
	assert( numCurves > 0 );

	double totLen = from.GetTotalLength();
	double lenPerSeg = totLen / numCurves;
	double timePerSeg = 1.f / double(numCurves);
	
	// make curves from knots :
	m_pieces.clear();
	m_pieces.resize(numCurves);
	
	Vec P[2];
	Vec D[2];
	int index=1;
	P[0] = from.GetValue(0.f);
	D[0] = from.GetDerivativeGlobal(0.f);

	for(int i = 0; i < numCurves; i++)
	{
		if ( i == numCurves-1 )
		{
			P[index] = from.GetValue(1.f);
			D[index] = from.GetDerivativeGlobal(1.f);
		}
		else
		{
			double accumLen = (i+1) * lenPerSeg;
		
			PlaceOnSpline place;
			from.GetPlaceByLengthAccurate(accumLen,&place);

			P[index] = from.GetValue(place);
			D[index] = from.GetDerivativeGlobal(place);
		}
		
		const Vec & P0 = P[index^1];
		const Vec & P1 = P[index];

		const Vec D0 = D[index^1] * timePerSeg;
		const Vec D1 = D[index] * timePerSeg;

		m_pieces[i].curve.SetFromEnds(P0,D0, P1,D1);
		m_pieces[i].accumTime = (i+1) * timePerSeg;

		index ^=1;
	}

	m_pieces[numCurves-1].accumTime = 1.f;
	
	m_type = eUniformLength;

	UpdateDerived();
}


void BezierSpline::SetSubdivisionSampled_RecursiveAdd( const Bezier & curve, const double t0, const double t1, const double maxErrSqr)
{
	if ( curve.ComputeLinearErrorBoundSqr() <= maxErrSqr )
	{
		 // ok, add it
		assert( m_pieces.size() > 0 || t0 == 0.f );
		assert( m_pieces.size() == 0 || m_pieces.back().accumTime == t0 );
		m_pieces.push_back(Piece());
		m_pieces.back().curve = curve;
		m_pieces.back().accumTime = t1;
	}
	else
	{
		// must subdivide
		Bezier b0,b1;
		curve.Subdivide(&b0,&b1);
		double tmid = (t0 + t1)*0.5f;

		SetSubdivisionSampled_RecursiveAdd(b0, t0,tmid, maxErrSqr );
		SetSubdivisionSampled_RecursiveAdd(b1, tmid,t1, maxErrSqr );
	}
}

void BezierSpline::SetSubdivisionSampled(const BezierSpline & from,const double maxErrSqr)
{
	m_pieces.clear();
	
	for(int i=0;i<from.GetNumSegments();i++)
	{
		double start,end;
		from.GetSegmentTimeRange(i,&start,&end);
		SetSubdivisionSampled_RecursiveAdd( from.m_pieces[i].curve, start,end, maxErrSqr);
	}
	
	m_type = eGeneric;

	UpdateDerived();
}

struct PieceToSplit
{
	double midTime; // refer to a piece by its midTime, since indexing is changing
	double errSqr;
	
	bool operator < (const PieceToSplit & rhs) const
	{
		return errSqr < rhs.errSqr;
	}
};

static double NumericalErrSqr(const BezierSpline &s1,const BezierSpline & s2,double start,double end,int numSteps)
{
	double errSqr = 0;
	
	// for two bezier curves we could integrate the exact error
	// but on the two splines the sampling might not be the same so we have to deal with various overlapping intervals
	//	definitely possible, but I'm lazy so just do some step samples :

	for(int i=0;i<numSteps;i++)
	{
		// no need to sample right at the ends, so offset in :
		double t = flerp(start,end,(double)(i+1)/(numSteps+1));
	
		Vec v1 = s1.GetValue(t);
		Vec v2 = s2.GetValue(t);
		
		errSqr += DistanceSqr(v1,v2);
	}
	
	// return err per step :
	errSqr /= numSteps;
	
	return errSqr;
}

static void InitPieceToSplit(PieceToSplit * pts,const BezierSpline & me,const BezierSpline & other,int piece)
{
	double start,end;
	me.GetSegmentTimeRange(piece,&start,&end);
	
	pts->midTime = faverage(start,end);
	
	// take a # of steps proportional to the time range
	const int c_totalNumSteps = 256;
	int numSteps = ftoi( c_totalNumSteps * (end - start) );
	numSteps = gClamped(numSteps,16,256);
	
	double errSqr = NumericalErrSqr(me,other,start,end,numSteps);
	
	// scale by time range so it's like the area of error
	//	@@ could scale by length
	errSqr *= (end - start);
	
	pts->errSqr = (double) errSqr;
}

void BezierSpline::MakeSimplification(const BezierSpline & from,const double maxErrSqr)
{
	// set from ends to start us with 1 segment
	SetUniformTimeSampled(from,2);
	
	// now keep splitting the worst one of my pieces until we're under error :
	
	Vector<PieceToSplit> heap;
	heap.push_back(PieceToSplit());
	InitPieceToSplit( &heap.back(), *this, from, 0 );
	
	while ( ! heap.empty() )
	{
		std::pop_heap(heap.begin(),heap.end());
		PieceToSplit cur = heap.back();
		heap.pop_back();
		
		if ( cur.errSqr < maxErrSqr )
			break; // done
			
		// find the piece we want to split :
		PlaceOnSpline place;
		GetPlaceByTime(cur.midTime,&place);
	
		// split place.segment
		int piece = place.segment;
		Piece toSplit = m_pieces[piece];
		
		double start,end;
		GetSegmentTimeRange(piece,&start,&end);
		
		m_pieces.insert( m_pieces.begin()+piece, toSplit );
		
		double tmid = faverage(start,end);
		
		Vec V0 = toSplit.curve.GetValue0();
		Vec D0 = toSplit.curve.GetDerivative0() * (1.f / (end - start));
		Vec V1 = toSplit.curve.GetValue1();
		Vec D1 = toSplit.curve.GetDerivative1() * (1.f / (end - start));
		
		/*assert( const Vec T0 = from.GetValue(start) );
		assert( const Vec T1 = from.GetValue(end) );
		assert( T0 == V0 );
		assert( T1 == V1 );*/
		
		// just add the new point at the middle
		// @@ in theory we could find the best spot to add a new point instead of just always doing the middle
		//  even just trying 1/3,1/2,2/3 might help a lot
		
		Vec VM = from.GetValue(tmid);
		Vec DM = from.GetDerivativeGlobal(tmid);
		
		// all derivatives are global, make them local :
		
		D0 *= (end - start)*0.5f;
		D1 *= (end - start)*0.5f;
		DM *= (end - start)*0.5f;
		
		m_pieces[piece].curve.SetFromEnds(V0,D0,VM,DM);
		m_pieces[piece].accumTime = tmid;
		
		m_pieces[piece+1].curve.SetFromEnds(VM,DM,V1,D1);
		m_pieces[piece+1].accumTime = end;
		
		heap.push_back(PieceToSplit());
		InitPieceToSplit( &heap.back(), *this, from, piece );
		std::push_heap(heap.begin(),heap.end());
		
		heap.push_back(PieceToSplit());
		InitPieceToSplit( &heap.back(), *this, from, piece+1 );
		std::push_heap(heap.begin(),heap.end());
	}
	
	UpdateDerived();
}

//! UpdateDerived uses only "curve"
void BezierSpline::UpdateDerived()
{
	if ( m_pieces.empty() )
		return;

	//m_bbox.SetToPoint( m_pieces[0].curve.GetValue0() );

	double accumLength = 0.f;

        for(int i = 0; i < (int)m_pieces.size(); i++)
	{
		double length = m_pieces[i].curve.ComputeLength();
		accumLength += length;
		m_pieces[i].accumLength = accumLength;
	}
}

//! ParameterizeByArcLength requires UpdateDerived
void BezierSpline::ParameterizeByArcLength()
{
	if ( m_pieces.empty() )
		return;

	double scale = 1.f / GetTotalLength();
        for(int i = 0; i< (int)m_pieces.size(); i++)
	{
		m_pieces[i].accumTime = m_pieces[i].accumLength * scale;
	}
	// make the ending time exactly 1.f
	assert( fisone(m_pieces.back().accumTime) );
	m_pieces.back().accumTime = 1.f;
}

double BezierSpline::GetTotalLength() const
{
	if ( m_pieces.empty() )
		return 0.f;
	return m_pieces.back().accumLength;
}

void BezierSpline::GetSegmentTimeRange(int seg,double *pStart,double *pEnd) const
{

        assert( seg >= 0 && seg < (int)m_pieces.size() );
	if ( seg == 0 )
		*pStart = 0.f;
	else
		*pStart = m_pieces[seg-1].accumTime;
	*pEnd = m_pieces[seg].accumTime;
}

void BezierSpline::GetSegmentLengthRange(int seg,double *pStart,double *pEnd) const
{
        assert( seg >= 0 && seg < (int)m_pieces.size() );
	if ( seg == 0 )
		*pStart = 0.f;
	else
		*pStart = m_pieces[seg-1].accumLength;
	*pEnd = m_pieces[seg].accumLength;
}

//! the derivative wrst global time is much larger
const Vec BezierSpline::GetDerivativeGlobal(const PlaceOnSpline & place) const
{
	double start,end;
	GetSegmentTimeRange(place.segment,&start,&end);
	double range = end-start;
	assert( range <= 1.f );
	Vec localD = GetCurve(place.segment).GetDerivative(place.localTime);
	return localD * (1.f/range);
}

void BezierSpline::GetPlaceByTime(const double t,PlaceOnSpline * pInto) const
{
	assert( fisinrange(t,0.f,1.f) );
	
	if ( m_pieces.empty() )
	{
		pInto->segment = 0;
		pInto->localTime = 0;
		return;
	}

	if( t <= 0.f )
	{
		pInto->segment = 0;
		pInto->localTime = 0;
		return;
	}

	if( t >= 1.f )
	{
		pInto->segment = m_pieces.size()-1;
		pInto->localTime = 1.f;
		return;
	}

	assert( m_pieces.back().accumTime == 1.f );

	if ( m_type == eUniformTime )
	{
		double v = t * m_pieces.size();
		pInto->segment = ftoi(v);
		pInto->localTime = v - pInto->segment;
	}
	else
	{
		// @@@@ TODO : could binary search on t for more speed
		int seg = 0;
                while ( t > m_pieces[seg].accumTime && seg < (int)m_pieces.size()-1 )
		{
			seg++;
		}

		// I'm in "seg"
		double start,end;
		GetSegmentTimeRange(seg,&start,&end);
		assert( fisinrange(t,start,end) );
		pInto->segment = seg;
		pInto->localTime = fmakelerpernoclamp(start,end,t);
	}
}

void BezierSpline::GetPlaceByLength(const double _len,PlaceOnSpline * pInto) const
{	
	if ( m_pieces.empty() )
	{
		pInto->segment = 0;
		pInto->localTime = 0;
		return;
	}

	if( _len <= 0.f )
	{
		pInto->segment = 0;
		pInto->localTime = 0;
		return;
	}

	if( _len >= GetTotalLength() )
	{
		pInto->segment = m_pieces.size()-1;
		pInto->localTime = 1.f;
		return;
	}

	const double len = fclamp(_len,0.f,GetTotalLength());
	
	if ( m_type == eUniformLength )
	{
		double v = len * m_pieces.size() / GetTotalLength();
		pInto->segment = ftoi(v);
		
		double start,end;
		GetSegmentLengthRange(pInto->segment,&start,&end);
		assert( fisinrange(len,start,end) );

		pInto->localTime = fmakelerpernoclamp(start,end,len);
	}
	else
	{
		// @@@@ TODO : could binary search on len for more speed
		int seg = 0;
		while ( len > m_pieces[seg].accumLength )
		{
			seg++;
		}

		// I'm in "seg"
		double start,end;
		GetSegmentLengthRange(seg,&start,&end);
		assert( fisinrange(len,start,end) );
		pInto->segment = seg;
		pInto->localTime = fmakelerpernoclamp(start,end,len);
	}
}

void BezierSpline::GetPlaceByLengthAccurate(const double _length,PlaceOnSpline * pInto) const
{
	const double length = fclamp(_length,0.f,GetTotalLength());
	// first get the seg :
	GetPlaceByLength(length,pInto);

	// I'm in "seg"
	double start,end;
	GetSegmentLengthRange(pInto->segment,&start,&end);
	assert( fisinrange(length,start,end) );

	double lenOnSeg = length - start;
	pInto->localTime = m_pieces[pInto->segment].curve.ComputeParameterForLength(lenOnSeg);
}

double BezierSpline::GetTimeFromPlace(const PlaceOnSpline & place) const
{
	if ( m_pieces.empty() )
		return 0.f;

	double start,end;
	GetSegmentTimeRange(place.segment,&start,&end);
	
	return flerp(start,end,place.localTime);
}

double BezierSpline::GetLengthFromPlace(const PlaceOnSpline & place) const
{
	if ( m_pieces.empty() )
		return 0.f;
	
	double start,end;
	GetSegmentLengthRange(place.segment,&start,&end);
	
	return flerp(start,end,place.localTime);
}

//! compare slines and return an error in distance squared
/*static*/ double BezierSpline::DifferenceSqr(const BezierSpline & s1,const BezierSpline & s2,const int numSamples)
{
	assert( numSamples > 2 );

	double worstErrSqr = 0.f;

	double dt = 1.f/(numSamples-1);
	for(int i=0;i<numSamples;i++)
	{
		double t = i * dt;
		const Vec v1 = s1.GetValue(t);
		const Vec v2 = s2.GetValue(t);
		double dSqr = DistanceSqr(v1,v2);
		worstErrSqr = gMAX(worstErrSqr,dSqr);
	}

	return worstErrSqr;
}

//! resample starting with "from" , to make a spline with as few controls as possible
//!	currently uses uniform-time sampling; eg. does not put fewer samples were the
//!	spline is smoother
void BezierSpline::MakeMinimumUniformTimeSampled(const BezierSpline & from,const double maxErrSqr)
{
	static const int c_numCompareSamples = 100;

	// do a binary search
	// will add 2 to searcher to make numKnots
	int searcher = from.GetNumSegments() -1;
	int searcherStep = (searcher+1)>>1;

	for( ;searcherStep > 1; searcherStep = (searcherStep+1)>>1)
	{
		int trySearcher = searcher - searcherStep;
		if ( trySearcher < 0 )
			continue;
		
		int tryKnots = trySearcher + 2;

		// @@ which sampler to use ?
		SetUniformTimeSampled(from,tryKnots);
		double errSqr = DifferenceSqr(*this,from,c_numCompareSamples);
		if ( errSqr <= maxErrSqr )
		{
			assert( trySearcher < searcher );
			searcher = trySearcher;
		}
	}

	int bestNumKnots = searcher+2;

	SetUniformTimeSampled(from,bestNumKnots);
}

bool BezierSpline::IsValid() const
{
	assert(m_type == eGeneric ||  m_type == eUniformTime ||  m_type == eUniformLength);
	return true;
}

double BezierSpline::GetTimeClosestPoint( const Vec & point ) const
{
	int n = m_pieces.size();

        double minDist = DBL_MAX;
	int minIndex = -1;

	// Find closest curve segment
        for(int i = 0; i < (int)m_pieces.size(); i++)
	{
		double dist = m_pieces[i].curve.MinDistanceToControlPoints(point);

		if(dist < minDist)
		{
			minDist = dist;
			minIndex = i;
		}
	}

	// Find time on closest segment
	double startTime, endTime;

	if(minIndex == 0)			
	{
		startTime = 0; 
		endTime = m_pieces[0].accumTime; 
	}
	else if (minIndex == n - 1)
	{ 
		startTime = m_pieces[Max(0, n - 2)].accumTime; 
		endTime = 1.0; 
	}
	else
	{
		startTime = m_pieces[minIndex - 1].accumTime;
		endTime = m_pieces[minIndex].accumTime;
	}

	// Find even closer point in current time period
	double closeTime = startTime;

	// Resolution
	double timeStep = (endTime - startTime) / 100;

	minDist = (point - GetValue(startTime)).norm();

	for(double t = startTime; t <= endTime; t += timeStep)
	{
		double dist = (point - GetValue(t)).norm();

		if(dist < minDist)
		{
			minDist = dist;
			closeTime = t;
		}
	}

	return closeTime;
}

Vec BezierSpline::ClosestPoint( const Vec & point ) const
{
	double t = GetTimeClosestPoint( point );

	return GetValue(t);
}

int IndexOf(const int i, const int size, const EEndHandling eEnds)
{
        if ( i < 0 )
        {
                if ( eEnds == eEnd_Clamp )
                        return 0;
                else // Wrap
                        return size + i;
        }
        else if ( i >= size )
        {
                if ( eEnds == eEnd_Clamp )
                        return size-1;
                else // Wrap
                        return i - size;
        }
        else
        {
                return i;
        }
}

