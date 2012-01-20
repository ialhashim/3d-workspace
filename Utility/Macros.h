#pragma once

// Standard C++
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "limits.h"
#include "float.h"
#include <time.h>
#include <string>
#include <numeric>

// GL extensions
#ifdef _WIN32
	#include <GL/GLee.h>
#else
	#include <GL/glew.h>
	#define GLEE_ARB_vertex_buffer_object GLEW_ARB_vertex_buffer_object
#endif

// Constants
static const double Epsilon_HIGH = DBL_EPSILON;
static const double Epsilon_LOW = FLT_EPSILON;
extern double Epsilon;

// Numerical stuff
#define FLOAT_INFINITY std::numeric_limits<float>::infinity()
#define DOUBLE_INFINITY std::numeric_limits<double>::infinity()

// STL containers
#include <vector>
#include <list>
#include <map>
#include <set>
#include <stack>

// STL beautification
#define StdVector std::vector
#define StdList std::list
#define StdMap std::map
#define StdSet std::set
#define Stack std::stack
#define StdString std::string
#define PairInt std::pair<int,int>
#define Pair std::pair
typedef StdVector < StdVector< float > > Vector2Df;

// IO name decoration
#define FileStream std::ifstream 
#define FileStreamOutput std::ofstream 
#define GetLine std::getline 

// Utility Macros
#define Max(a,b) (((a) > (b)) ? (a) : (b))
#define Min(a,b) (((a) < (b)) ? (a) : (b))
#define MaxOf(a,b,c) (Max(a, Max(b,c)))
#define Mod(x,m) ((x % m + m) % m)
#define RANGED(min, v, max) ( Max(min, Min(v, max)) ) 
#define RANGE(i, min, max) (  ((i >= min) && (i <= max)) ? 1 : 0)
#define SIGN(i) ((i >= 0) ? (1) : (-1))
#define SET_HAS(SET, OBJECT) ((SET.find(OBJECT) != SET.end()) ? 1 : 0)
#define MAP_HAS(MAP, OBJECT) ((MAP.find(OBJECT) != MAP.end()) ? 1 : 0)
#define PREV(i, N) ((i + N-1) % N)
#define NEXT(i, N) ((i + 1) % N)
#define SWAP(x, y, T) do { T temp##x##y = x; x = y; y = temp##x##y; } while (0)
#define AROUND(x, target, threshold) ( (abs(x) - abs(target) < threshold) ? 1 : 0)
#define RADIANS(deg)    ((deg)/180.0 * M_PI)
#define DEGREES(rad)    ((rad)/M_PI * 180.0)
#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

// Basic STL converters
template <typename T> static inline StdVector<T> SET_TO_VECTOR(StdSet<T> fromSet){
	StdVector<T> result;
        for(typename StdSet<T>::iterator it = fromSet.begin(); it != fromSet.end(); it++)
		result.push_back(*it);
	return result;
}
template <typename T> static inline StdSet<T> VECTOR_TO_SET(StdVector<T> fromVector){
	StdSet<T> result;
        for(typename StdVector<T>::iterator it = fromVector.begin(); it != fromVector.end(); it++)
		result.insert(*it);
	return result;
}
template <typename T> static inline StdList<T> VECTOR_TO_LIST(StdVector<T> fromVector){
	StdList<T> result;
        for(typename StdVector<T>::iterator it = fromVector.begin(); it != fromVector.end(); it++)
		result.push_back(*it);
	return result;
}
template <typename T> static inline StdVector<T> LIST_TO_VECTOR(StdList<T> fromList){
	StdVector<T> result;
        for(typename StdList<T>::iterator it = fromList.begin(); it != fromList.end(); it++)
		result.push_back(*it);
	return result;
}
template <typename T> static inline StdSet<T> LIST_TO_SET(StdList<T> fromList){
	StdSet<T> result;
        for(typename StdList<T>::iterator it = fromList.begin(); it != fromList.end(); it++)
		result.insert(*it);
	return result;
}

double inline uniform(double a = 0.0, double b = 1.0)
{
	double len = b - a;
	return ((double)rand()/RAND_MAX) * len + a;
}

unsigned inline int fact(unsigned int n){
	unsigned int i, p=1;

	for(i = 2; i <= n; i++) 
		p *= i;

	return p;
}

double inline gaussianFunction(double x, double mu = 0.0, double sigma = 1.0){
	//double a = 1.0 / (sigma * sqrt(2 * M_PI));
	double b = mu;
	double c = sigma;

	// normalized guassian with N(0, \sigma) = 1
	return exp( - (pow(x - b, 2) / (2 * pow(c, 2)) ) );
}

// Rodrigues' rotation
#define ROTATE_VEC(v, theta, axis) (v = v * cos(theta) + cross(axis, v) * sin(theta) + axis * dot(axis, v) * (1 - cos(theta)))
#define ROTATED_VEC(v, theta, axis) (v * cos(theta) + cross(axis, v) * sin(theta) + axis * dot(axis, v) * (1 - cos(theta)))

template <typename VECTYPE>
VECTYPE inline RotateAround(VECTYPE point, VECTYPE pivot, VECTYPE axis, double theta){
	VECTYPE result = point;

	result -= pivot;
	result = ROTATE_VEC(result, theta, axis);
	result += pivot;

	return result;
}

template <typename VECTYPE>
void inline RotateFromTo(VECTYPE from, VECTYPE to, VECTYPE & point, VECTYPE pivot = VECTYPE(0,0,0))
{
	VECTYPE axis = cross(from, to).normalized();
	double theta = acos(dot(from.normalize(), to.normalize()));

	point -= pivot;
	point = ROTATE_VEC(point, theta, axis);
	point += pivot;
}

// Array operations
#include <algorithm>
#define MaxElement(v) (*max_element(v.begin(), v.end()))
#define MinElement(v) (*min_element(v.begin(), v.end()))
#define DivideVector(a,value) for(int i = 0; i < a.size(); ++i) a[i] /= value;
#define Sum(v)  (std::accumulate(v.begin(), v.end(), 0.0))
#define Avg(v)  (Sum(v) / v.size())

// Timer
#include <QElapsedTimer>
#define Timer QElapsedTimer
#define CreateTimer(timer)  QElapsedTimer timer; timer.start()
