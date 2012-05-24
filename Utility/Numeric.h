#pragma once

#include <vector>
#include "GraphicsLibrary/Mesh/SurfaceMesh/Vector.h"
#include <QColor>
#include <QString>

#include <Eigen/Dense>

extern double GC_GAUSSIAN_SIGMA;

typedef std::vector< std::vector<double> >	Buffer2d;
typedef std::vector< std::vector<bool> >	Buffer2b;
typedef std::vector< std::vector<Vec2i> >   Buffer2v2i;

typedef Vec3d Point;

// Extrema
double getMinValue( Buffer2d & image );
double getMaxValue( Buffer2d & image );
Vec3d maxMidMinValues( Buffer2d & image );

// Region
double maxValueInRegion( Buffer2d& image,  std::vector< Vec2i >& region);
Vec2i sizeofRegion( std::vector< Vec2i >& region );
void BBofRegion( std::vector< Vec2i >& region, Vec2i &bbmin, Vec2i &bbmax );
Vec2i centerOfRegion( std::vector< Vec2i >& region );
std::vector< double > getValuesInRegion( Buffer2d& image, std::vector< Vec2i >& region, bool xFlipped = false );
std::vector< Vec2i > getRegionGreaterThan( Buffer2d& image, Buffer2b& mask, Vec2i seed, double threshold);
std::vector< std::vector< Vec2i > > getRegionsGreaterThan(Buffer2d& image, double threshold);

// Shifting
std::vector< Vec2i > deltaVectorsToKRing(int deltaX, int deltaY, int K);
std::vector< Vec2i > shiftRegionInBB( std::vector< Vec2i >& region, Vec2i delta, Vec2i bbmin, Vec2i bbmax );

// Color
QRgb jetColor( double val, double min, double max );
void setRegionColor( Buffer2d& image, std::vector< Vec2i >& region, double color );
void setPixelColor( Buffer2d& image, Vec2i pos, double color );

// Visualization
void visualizeRegions( int w, int h, std::vector< std::vector<Vec2i> >& regions, QString filename );
template< typename T >
std::vector< std::vector < T > > createImage( int w, int h, T intial);
void saveAsBinaryImage( Buffer2d& image, QString fileName );
void saveAsImage( Buffer2d& image, QString fileName );
void saveAsData( Buffer2d& image, double maxV, QString fileName );

// Random set of unique integers
std::vector<int> uniformIntegerSamples(int N, int range);
std::vector<Vec2i> sampleRegion(std::vector<Vec2i> &region, int N);

// Rotation
Eigen::Matrix3d rotationMatrixAroundAxis(Vec3d u, double theta);
Vec3d rotatePointByMatrix( Eigen::Matrix3d &R, Vec3d p );

// Cuboid volume
double volumeOfBB(Vec3d &extents);

// Points
void twoFurthestPoints( std::vector<Point> &points, Point &p1, Point &p2 );
double distanceCluster2Cluster( std::vector<Point> &cluster1, std::vector<Point> &cluster2 );
std::vector< std::vector<Point> > twoClassClustering( std::vector<Point>& points, Point seed1, Point seed2 );
Point centroid(std::vector<Point> points);

// Sample poly-line
std::vector<Point> uniformSampleCurve(std::vector<Point> & points);

// Adaptive maximum region detection
Buffer2v2i getMaximumRegions(Buffer2d &image);

// AABB
std::vector<Point> cornersOfAABB(Vec3d bbmin, Vec3d bbmax);
void computeAABB( std::vector<Point> &points, Vec3d &bbmin, Vec3d &bbmax);