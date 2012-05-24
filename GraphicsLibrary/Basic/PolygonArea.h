#pragma once
// Functions from paper entitled "Fast Polygon Area and Newell Normal Computation",
// Daniel Sunday, journal of graphics tools, 7(2):9-13, 2002.
// http://www.acm.org/jgt/papers/Sunday02/FastArea.html

// Input: 2D polygon
static inline double findArea(int n, double *x, double *y)        
{
	// Assume that the 2D polygon's vertex coordinates are
	// stored in (x[0], y[0]), ..., (x[n-1], y[n-1]),
	// with no assumed vertex replication at the end.

	// Initialize sum with the boundary vertex terms
	double sum = x[0] * (y[1] - y[n-1]) + x[n-1] * (y[0] - y[n-2]);

	for (int i=1; i < n-1; i++) {
		sum += x[i] * ( y[i+1] - y[i-1] );
	}

	return (sum / 2.0);
}

// return the signed area of a 3D planar polygon (given normal vector)
// 3D planar polygon, and plane normal
static inline double findArea3D(int n, double *x, double *y, double *z, double nx, double ny, double nz) 
{
	// select largest normal coordinate to ignore for projection
	double ax = (nx>0 ? nx : -nx);	// abs nx
	double ay = (ny>0 ? ny : -ny);	// abs ny
	double az = (nz>0 ? nz : -nz);	// abs nz
	double len = sqrt(nx*nx + ny*ny + nz*nz); // length of normal

	if (ax > ay) {
		if (ax > az)			       // ignore x-coord
			return findArea(n, y, z) * (len / nx);
	}
	else if (ay > az)			       // ignore y-coord
		return findArea(n, z, x) * (len / ny);
	return findArea(n, x, y) * (len / nz); // ignore z-coord
}

// output the approximate unit normal of a 3D nearly planar polygon
// return the area of the polygon
// Inputs: 3D polygon, output unit normal
static inline double findNormal3D(int n, double *x, double *y, double *z, double *nx, double *ny, double *nz) 
{
	// get the Newell normal
	double nwx = findArea(n, y, z);
	double nwy = findArea(n, z, x);
	double nwz = findArea(n, x, y);
	// get length of the Newell normal
	double nlen = sqrt( nwx*nwx + nwy*nwy + nwz*nwz );
	// compute the unit normal
	*nx = nwx / nlen;
	*ny = nwy / nlen;
	*nz = nwz / nlen;
	return nlen;    // area of polygon = length of Newell normal
}

static inline Vec3d centerOfPoints(const std::vector<Vec3d>& points)
{
	Vec3d center;

	for(std::vector<Vec3d>::const_iterator it = points.begin(); it != points.end(); it++)
		center += *it;

	center /= points.size();

	return center;
}

static inline double findArea3D(const std::vector<Vec3d>& points, const Vec3d& n)
{
	int N = points.size();

	std::vector<double> x(N), y(N), z(N);

	for(int i = 0; i < N; i++)
	{
		x[i] = points[i].x();
		y[i] = points[i].y();
		z[i] = points[i].z();
	}

	double * X = &x.front();
	double * Y = &y.front();
	double * Z = &z.front();

	double area = findArea3D(N, X, Y, Z, n.x(), n.y(), n.z());

	return area;
}

static inline double findNormal3D(const std::vector<Vec3d>& points, Vec3d& n)
{
	int N = points.size();

	std::vector<double> x(N), y(N), z(N);

	for(int i = 0; i < N; i++)
	{
		x[i] = points[i].x();
		y[i] = points[i].y();
		z[i] = points[i].z();
	}

	double * X = &x.front();
	double * Y = &y.front();
	double * Z = &z.front();

	return findNormal3D(N, X, Y, Z, &n.x(), &n.y(), &n.z());
}

static inline double signedArea(const std::vector<Vec3d>& points, const Vec3d& n, const Vec3d& center)
{
	int N = points.size();
	Vec3d sumVec;

	for(int i = 0; i < N; i++)
	{
		int j = NEXT(i, N);
		sumVec += cross((points[i]-center),(points[j]-center));
	}

	return dot(sumVec, n);
}
