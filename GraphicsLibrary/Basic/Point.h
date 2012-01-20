#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Utility/Macros.h"

#include "GUI/Viewer/libQGLViewer/QGLViewer/qglviewer.h"
using namespace qglviewer;

struct Point2D
{
public:
	int x;
	int y;

	Point2D(){x = y = 0;}
	Point2D(int x, int y)
	{
		this->x = x;
		this->y = y;
	}
	Point2D(const Point2D& fromPoint)
	{
		this->x = fromPoint.x;
		this->y = fromPoint.y;
	}

	friend bool operator< (const Point2D &a, const Point2D &b)
	{
		return (a.x < b.x) || ( (a.x == b.x) && (a.y < b.y) );
	}

	Point2D operator* (const int scale)
	{
		Point2D result(*this);
		result.x *= scale;
		result.y *= scale;
		return result;
	}

	Point2D operator+ (const Point2D delta)
	{
		Point2D result(*this);
		result.x += delta.x;
		result.y += delta.y;
		return result;
	}

	Point2D& operator= (const Point2D& fromPoint) 
	{
		this->x = fromPoint.x;
		this->y = fromPoint.y;
		return *this;
	}

	inline bool operator== (const Point2D& p) const
	{
		return (this->x == p.x && this->y == p.y);
	}

	inline void hash_combine(size_t& seed, int v) const
	{
		seed ^= static_cast<size_t>(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
	}

	operator size_t() const
	{
		size_t seed = 0;
		hash_combine(seed, x);
		hash_combine(seed, y);
		return seed;
	}

	inline void add(int deltaX, int deltaY)
	{
		x += deltaX;
		y += deltaY;
	}

	inline double distanceTo(const Point2D& otherPoint)
	{
		return sqrt( pow(2.0, x - otherPoint.x) + pow(2.0, y - otherPoint.y) );
	}

	inline double DistSq() const { return x*x + y*y; }
};

#define Vector2DPoint std::StdVector<std::StdVector<Point2D> >

class Point3D : public Vec{

public:

	Point3D() : Vec() {}

	Point3D(double fromX, double fromY, double fromZ){
		x = fromX;
		y = fromY;
		z = fromZ;
	}

	Point3D(const Vec& from){
		x = from.x;
		y = from.y;
		z = from.z;
	}

	void set(const double& X, const double& Y, const double& Z){
		x = X;
		y = Y;
		z = Z;
	}

	void set(const Vec& from){
		x = from.x;
		y = from.y;
		z = from.z;
	}

	void set(const Point3D& from){
		x = from.x;
		y = from.y;
		z = from.z;
	}

	void add(const double& X, const double& Y, const double& Z){
		x += X;
		y += Y;
		z += Z;
	}

	int above(Vec & pos, Vec & dir){
		double sign = dir * (*this - pos);

		if(sign > 0)
			return 1;
		else if(sign < 0)
			return -1;
		else
			return 0;
	}

	static double cotangent(Point3D & v0, Point3D & v1){
		double f = (v0 ^ v1).norm();
		if(f == 0) return 0;
		return (v0 * v1) / f;
	}

	const Vec& vec(){
		return *this;
	}

	static Vec interpolateVec(double t, const Vec & from, const Vec & to)
	{
		if(t == 0.0)	return from;
		if(t == 1.0)	return to;

		return from + (t * (to - from));
	}

	static Vec MidVec(const Vec & a, const Vec & b)
	{
		return Vec((a.x + b.x) / 2.0, (a.y + b.y) / 2.0, (a.z + b.z) / 2.0);
	}

	static Vec MidVec(const StdVector<Vec> & points)
	{
		int N = points.size();

		Vec sum;

		for(StdVector<Vec>::const_iterator it = points.begin(); it != points.end(); it++)
			sum += *it;

		return sum / (double)N;
	}

	static double maxLength(Point3D &a, Point3D &b, Point3D &c)
	{
		return Max((a-b).norm(), Max((a-c).norm(), (b-c).norm()));
	}

	static double minLength(Point3D &a, Point3D &b, Point3D &c)
	{
		return Min((a-b).norm(), Min((a-c).norm(), (b-c).norm()));
	}

	static double totalLength(Point3D &a, Point3D &b, Point3D &c)
	{
		return (a-b).norm() + (a-c).norm() + (b-c).norm();
	}

	static Vec RotateAround(const Vec & point, const Vec & pivot, const qglviewer::Quaternion & q)
	{
		Vec result = point;

		result -= pivot;
		result = q.rotate(result);
		result += pivot;

		return result;
	}

	static double angleBetweenTwo(const Vec &a, const Vec &b)
	{
		double Dot = (a).unit() * (b).unit();

		if(abs(Dot) >= 1.0)	
			Dot = 1.0; 

		return acos(Dot);
	}

	double angleBetween(const Vec &a, const Vec &b)
	{
		double Dot = (a - *this).unit() * (b - *this).unit();

		if(abs(Dot) >= 1.0)	
			Dot = 1.0; 

		return acos(Dot);
	}

	static double minAngle(const Vec &a, const Vec &b, const Vec &c)
	{
		return Min(Point3D(a).angleBetween(b,c), 
			Min(Point3D(b).angleBetween(c,a), 
			Point3D(c).angleBetween(a,b)));
	}

	static double maxAngle(const Vec &a, const Vec &b, const Vec &c)
	{
		return Max(Point3D(a).angleBetween(b,c), 
			Max(Point3D(b).angleBetween(c,a), 
			Point3D(c).angleBetween(a,b)));
	}

	static double halfAlphaTangent(const Vec &a, const Vec &b, const Vec &c)
	{
		double Dot = (a - b).unit() * (c - b).unit();

		Dot = RANGED(-1.0, Dot, 1.0); // bound

		return  tan(acos(Dot) / 2.0); // tan (angle / 2)
	}

	static double singed_angle(const Vec &a, const Vec &b, const Vec &axis)
	{
		double cosAngle = a.unit() * b.unit();
		double angle = acos( RANGED(-1.0, cosAngle, 1.0) );

		Vec c = a ^ b;

		if (c * axis < 0)
			return -angle;

		return angle;
	}

	static bool isBetween(const Vec& p, const Vec& a, const Vec& b, double Eps = 0.0)
	{
		double length = (a - b).norm();
		double dist = (p - a).norm() + (p - b).norm();

		if(dist > length + Eps)
			return false;
		else
			return true;
	}

	static bool isSameSide(const Vec& p1, const Vec& p2, const Vec& a, const Vec& b)
	{
		Vec cp1 = (b-a) ^ (p1-a);
		Vec cp2 = (b-a) ^ (p2-a);

		return (cp1 * cp2 >= 0);
	}

	StdString toString()
	{
		char buff[2048];

		sprintf(buff,"v %f %f %f", x, y, z);

		return buff;
	}

	void OrthogonalBasis( Vec &left, Vec &up ) const
	{
		float l, s;

		if ( fabs( z ) > 0.7f ) {
			l = y * y + z * z;
			s = 1.0f / sqrt ( l );
			up[0] = 0;
			up[1] = z * s;
			up[2] = -y * s;
			left[0] = l * s;
			left[1] = -x * up[2];
			left[2] = x * up[1];
		}
		else {
			l = x * x + y * y;
			s = 1.0f / sqrt ( l );
			left[0] = -y * s;
			left[1] = x * s;
			left[2] = 0;
			up[0] = -z * left[1];
			up[1] = z * left[0];
			up[2] = l * s;
		}
	}
};

struct VecIndex
{
	Vec v;
	int index;

	VecIndex(){	index = -1; }
	VecIndex(Vec vec, int Index) : v(vec), index(Index){}
	operator const Vec() const {	return v;	}
};

struct Rect
{
	int l,r;
	int t,b;

	Rect(int Left, int Right, int Top, int Bottom)
	{l = Left; r = Right; t = Top; b = Bottom;}
	Rect(){l = r = t = b = 0;}

	int width() const{return abs(r - l);}
	int height() const{return abs(b - t);}
};
