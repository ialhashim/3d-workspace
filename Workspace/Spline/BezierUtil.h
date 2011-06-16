#pragma once

#include <float.h> // for FLT_MAX
#include <math.h>
#include <malloc.h> 
#include <assert.h>
#include <algorithm>
#include <string.h> // for memcpy

#include "Point.h"

// Based on Charles Bloom Galaxy3 http://www.cbloom.com/3d/galaxy3/index.html

/*
The Bloom Public License (BPL)

Copyright (c) 1998-2009, Charles Bloom

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

* Neither the name Charles Bloom, "cbloom" nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

Any usage allowed by the GNU Public License (GPL) is allowed.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

//-------------------------------------------------------------------
// float consts :

#undef PI
#undef TWO_PI
#undef HALF_PI
#undef EPSILON

static const float PI                   = 3.141592653590f;
static const float TWO_PI               = 6.2831853071796f;
static const float HALF_PI              = 1.570796326795f;

static const float EPSILON              = 1e-3f; // @@ scary; much more thought needed

static const float RADS_PER_DEG = 0.01745329251994329576f;
static const float DEGS_PER_RAD = 57.2957795130823208767f;

//@@ get rid of these or the non-g versions
#define gPI             PI
#define gTWO_PI         TWO_PI
#define gHALF_PI        HALF_PI
#define gEPSILON        EPSILON

#define FLOAT_AS_INT(f)         (reinterpret_cast<const ulong &>(f))
#define FLOAT_AS_INT_MUTABLE(f) (reinterpret_cast<ulong &>(f))
#define INT_AS_FLOAT(i)         (reinterpret_cast<const float &>(i))
#define INT_AS_FLOAT_MUTABLE(i) (reinterpret_cast<float &>(i))

template<class Type> inline void gSwap(Type &a,Type &b)
{
	Type c = a; a = b; b = c;
}

template<class Type,class Type1,class Type2> inline Type gClamped(const Type &x,const Type1 &lo,const Type2 &hi)
{
	return ( ( x < lo ) ? lo : ( x > hi ) ? hi : x );
}

#define gMIN(a,b)       ( ((a) < (b)) ? (a) : (b) )
#define gMAX(a,b)       ( ((a) < (b)) ? (b) : (a) )

#define gMIN3(a,b,c)    gMIN(gMIN(a,b),c)
#define gMAX3(a,b,c)    gMAX(gMAX(a,b),c)

#define gABS(a)         ( ((a) < 0) ? -(a) : (a) )

//-------------------------------------------------------------------
// float utility functions :

inline bool fisvalid(const float f)
{
	// this works because NAN always returns false on compares :
	return (f >= -FLT_MAX && f <= FLT_MAX );
}

//! must do this lame shit to make sure both values are
//      treated the same; if not, one is left in a register and
//      is actually a double
bool fequalexact(const float x,const float y);

//! float == test with slop
inline bool fequal(const float f1,const float f2,const float tolerance = EPSILON)
{
	const float diff = fabsf(f1 - f2);
	return diff <= tolerance;
}

//! return a clamped float in the range lo to hi
inline float fclamp(const float x,const float lo,const float hi)
{
	return ( ( x < lo ) ? lo : ( x > hi ) ? hi : x );
}

//! return a clamped float in the range 0 to 1
inline float fclampunit(const float x)
{
	return ( ( x < 0.f ) ? 0.f : ( x > 1.f ) ? 1.f : x );
}

inline float faverage(const float x,const float y)
{
	return (x+y)*0.5f;
}

//! return a lerped float at time "t" in the interval from lo to hi
//      t need not be in [0,1] , we can extrapolate too
inline float flerp(const float lo,const float hi,const float t)
{
	return lo + t * (hi - lo);
}

inline float flerpangle(float fm, float to, float t)
{
	// make sure we take the short way around the circle
	while( fm > (to + PI) )
	{
		to += TWO_PI;
	}
	while( fm < (to - PI) )
	{
		to -= TWO_PI;
	}
	return flerp(fm,to,t);
}

//! return a float such that flerp(lo,hi, fmakelerper(lo,hi,v) ) == v
inline float fmakelerpernoclamp(const float lo,const float hi,const float v)
{
	assert( hi > lo );
	return (v - lo) / (hi - lo);
}

inline float fmakelerperclamp(const float lo,const float hi,const float v)
{
	assert( hi > lo );
	return fclampunit( (v - lo) / (hi - lo) );
}

//!     remap (0->1) to (0->1) but with derivative continuity
inline float fHermiteLerpParameter(const float t)
{
	// doesn't make sense outside of [0,1]
	assert( t >= - EPSILON && t <= 1.f + EPSILON );
	return (3.f - 2.f*t)*t*t;
}

//! use fHermiteLerpClamping if t could be out of the [0,1] range
inline float fHermiteLerpClamping(const float t)
{
	assert( fisvalid(t) );
	if ( t <= 0.f ) return 0.f;
	else if ( t >= 1.f ) return 1.f;
	else return fHermiteLerpParameter(t);
}

//! even better continuity than Hermite lerp, but slower
inline float fCosLerpParameter(const float t)
{
	// doesn't make sense outside of [0,1]
	assert( t >= - EPSILON && t <= 1.f + EPSILON );
	return 0.5f - 0.5f * cosf( t * PI );
}

//! use fCosLerpClamping if t could be out of the [0,1] range
inline float fCosLerpClamping(const float t)
{
	assert( fisvalid(t) );
	if ( t <= 0.f ) return 0.f;
	else if ( t >= 1.f ) return 1.f;
	else return fCosLerpParameter(t);
}

//! fsign; returns 1.f or -1.f for positive or negative float
//! could do a much faster version using FLOAT_AS_INT if this is needed
inline float fsign(const float f)
{
	return (f >= 0.f) ? 1.f : -1.f;
}

// handles 0's to err on the side of returning more true's
inline bool fsignsame(const float a,const float b)
{
	if ( a <= 0.f && b <= 0.f ) return true;
	if ( a >= 0.f && b >= 0.f ) return true;
	return false;
}

// clamp f to positive
inline float fpos(const float f)
{
	if ( f < 0.f ) return 0.f;
	return f;
}

// clamp f to negative
inline float fneg(const float f)
{
	if ( f > 0.f ) return 0.f;
	return f;
}

inline float fsquare(const float f)
{
	return f*f;
}

inline float fcube(const float f)
{
	return f*f*f;
}

//! float == 0.f test with slop
inline bool fiszero(const float f1,const float tolerance = EPSILON)
{
	return fabsf(f1) <= tolerance;
}

//! return (f1 in [0,1]) with slop
inline bool fiszerotoone(const float f1,const float tolerance = EPSILON)
{
	return (f1 >= -tolerance && f1 <= (1.f + tolerance) );
}

//inline bool fisinrange_exact(const float f,const float lo,const float hi)
inline bool fisinrange(const float f,const float lo,const float hi)
{
	return f >= lo && f <= hi;
}

inline bool fisinrange_eps(const float f,const float lo,const float hi,const float tolerance = EPSILON)
{
	return f >= lo-tolerance && f <= hi+tolerance;
}

//! float == 1.f test with slop
inline bool fisone(const float f1,const float tolerance = EPSILON)
{
	return fabsf(f1-1.f) <= tolerance;
}

//! float to int is expensive
//! this depends on the FPU rounding mode
//! asserts if the rounding mode is not C-rounding compatible
#if 1 //{
inline int ftoi(const float f)
{
        return (int) f;
}
inline int dtoi(const double f)
{
        return (int) f;
}
#else //}{
inline int ftoi(const float f)
{
	return int(f);
}
inline int dtoi(const double f)
{
	return int(f);
}
#endif //}

//! int to float is easy
inline float itof(const int i)
{
	return (float)i;
}

inline int intlog2(const int i)
{
	const float f = float(i);
	return intlog2(f);
}


//! degrees to radians converter
inline float fdegtorad(const float d)
{
	return d * RADS_PER_DEG;
}
//! radians to degrees converter
inline float fradtodeg(const float r)
{
	return r * DEGS_PER_RAD;
}

inline void finvalidate(float & f)
{
	FLOAT_AS_INT_MUTABLE(f) = 0xFFFFFFFF;
	//assert( ! fisvalid(f) );
}

/*
QuadraticSplineLerper

makes a spline P(t) which has Bezier controls
B0 = {0,0}
B1 = "point"
B2 = {1,1}

if "point" as {0.5,0.5} it would be a straight line lerp
We solve for

P_y(P_x)

"point" must be in range so that the solution is single-valued
*/
inline float fquadraticsplinelerper(const float f, const float ptX, const float ptY )
{
	assert( fisinrange(ptX,0.f,1.f) );
	assert( fisinrange(ptY,0.f,1.f) );
	assert( !fequal(ptX,0.5f) );    //singularity here.

	float bx = ptX;
	float a = (1 - 2.f*bx);
	float A = bx*bx + a*f;
	float t = (sqrtf(A) - bx)/a;
	float y = t*t + ptY*2.f*t*(1.f-t);
	return y;
}

//-------------------------------------------------------------------
// little utilities to warp the unit range 0->1 to itself

//! simple squaring warps x down
//!     result >= input (equality at 0 and 1)
inline float fwarpunit_down(const float x)
{
	return x*x;
}

//! cubic polynomial
//!     result <= input (equality at 0 and 1)
inline float fwarpunit_up(const float x)
{
	const float y = (1.f - x);

	// the 1.5 in here is actually 4 * ( 2 * k - 1 )
	//      where "k" is the value of the curve at x = 0.5
	// 1.5 comes from k = 0.69 , roughly, which is near sqrt(.5)
	// 1.65688 actually gives you sqrt(.5)

	return x * x * (x + 3.f * y) + 1.5f * x * y * y;
}

//-------------------------------------------------------------------

// min/max of absolute values. returns the original value
//
inline float fmaxabs( const float a, const float b )
{
	return fabsf(a) > fabsf(b) ? a : b;
}

inline float fminabs( const float a, const float b )
{
	return fabsf(a) < fabsf(b) ? a : b;
}

//-------------------------------------------------------------------

//! solves the equations of dynamics :
//!     (d/dt) val = (1 / time_scale) * (towards - val);
//!     that is, "val" goes towards "towards" in a time scale of "time_scale"
//!     this is linear drag with an external force
inline float DampedDrive(const float val, const float towards, const float time_scale, const float time_step)
{
	// amazingly enough this becomes just a non-linear lerp
	//      for (time_step / time_scale) small, it becomes a linear lerp !
	const float z = expf( - time_step / time_scale );
	return val * z + towards * (1.f - z);
}

inline float LerpedDrive(const float from, const float to, const float speed, const float time_step)
{
	const float step = speed * time_step;
	const float delta = to - from;
	if ( delta >= 0.f )
	{
		if ( delta <= step )
		{
			return to;
		}
		else
		{
			return from + step;
		}
	}
	else
	{
		if ( delta >= -step )
		{
			return to;
		}
		else
		{
			return from - step;
		}
	}
}

//-------------------------------------------------------------------
/***

(use of these is transparent to the user)

#define to replace some of the standard C math functions
with versions that will assert about the range of arguments

Make these #define *before* the safe functions

***/

#ifdef DO_assertS //{

extern float sqrtf_asserting(const float f);
extern float asinf_asserting(const float f);
extern float acosf_asserting(const float f);
extern double sqrt_asserting(const double f);
extern double asin_asserting(const double f);
extern double acos_asserting(const double f);

#define sqrtf   sqrtf_asserting
#define asinf   asinf_asserting
#define acosf   acosf_asserting
#define sqrt    sqrt_asserting
#define asin    asin_asserting
#define acos    acos_asserting

#endif //} DO_assertS

//-------------------------------------------------------------------
/***

These assert that the arguments are loosely in range, but not
perfectly.  The idea is for things like this :

v.Normalize();
angle = acosf(v.x);

So, v.x should be in (-1 to 1), but cuz of epsilons, it could be
slightly out of range.  So, typically you'll just start with the
above code.  Then, hopefully, the _asserting code will catch it
and fire.  Then you can fix it up by doing this :

v.Normalize();
angle = acosf_safe(v.x);

****/

inline float sqrtf_safe(const float f)
{
	// f should be >= 0.f , but allow a tiny negative & clamp it
	assert( f >= - EPSILON );
	if ( f <= 0.f )
		return 0.f;
	else
		return (float) sqrtf(f);
}

inline float asinf_safe(const float f)
{
	// f should be in [-1,1] but allow a little slip & clamp it
	assert( f >= (- 1.f - EPSILON) && f <= (1.f + EPSILON) );

	const float fclamped = fclamp(f,-1.f,1.f);

	return (float) asinf(fclamped);
}

inline float acosf_safe(const float f)
{
	// f should be in [-1,1] but allow a little slip & clamp it
	assert( f >= (- 1.f - EPSILON) && f <= (1.f + EPSILON) );

	const float fclamped = fclamp(f,-1.f,1.f);

	return (float) acosf(fclamped);
}

inline double sqrt_safe(const double f)
{
	// f should be >= 0.f , but allow a tiny negative & clamp it
	assert( f >= - EPSILON );
	if ( f <= 0. )
		return 0.;
	else
		return sqrt(f);
}

inline double asin_safe(const double f)
{
	// f should be in [-1,1] but allow a little slip & clamp it
	assert( f >= (- 1. - EPSILON) && f <= (1. + EPSILON) );

	const double fclamped = gClamped(f,-1.,1.);

	return asin(fclamped);
}

inline double acos_safe(const double f)
{
	// f should be in [-1,1] but allow a little slip & clamp it
	assert( f >= (- 1.f - EPSILON) && f <= (1.f + EPSILON) );

	const double fclamped = gClamped(f,-1.,1.);

	return acos(fclamped);
}

//-------------------------------------------------------------------
inline const Vec MakeLerp(const Vec &v1, const Vec &v2, const float t)
{
	return Vec(
		v1.x + t * (v2.x - v1.x),
		v1.y + t * (v2.y - v1.y),
		v1.z + t * (v2.z - v1.z)
		);
}

inline float DistanceSqr(const Vec & a, const Vec &b)
{
	return fsquare(a.x - b.x) + fsquare(a.y - b.y) + fsquare(a.z - b.z);
}


inline float Distance(const Vec & a,const Vec &b)
{
	return sqrtf( DistanceSqr(a,b) );
}

inline void SetAverage(Vec & result, const Vec &v1, const Vec &v2)
{
	// gWarnING : may produce different results than SetLerp(0.5) due to floating round-off
	result.x = 0.5f * (v2.x + v1.x);
	result.y = 0.5f * (v2.y + v1.y);
	result.z = 0.5f * (v2.z + v1.z);
}

inline void SetLerp(Vec & result, const Vec &v1, const Vec &v2, const float t)
{
	result.x = v1.x + t * (v2.x - v1.x);
	result.y = v1.y + t * (v2.y - v1.y);
	result.z = v1.z + t * (v2.z - v1.z);
}
//-------------------------------------------------------------------
