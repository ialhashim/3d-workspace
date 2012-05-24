#include "Numeric.h"

#include "Utility/Macros.h"
#include "Utility/ColorMap.h"
#include <QRect>
#include <QImage>
#include <QFile>
#include <stack>

double GC_GAUSSIAN_SIGMA = 0.2;


// Extreme
double getMaxValue( Buffer2d& image )
{
	int h = image.size();

	std::vector< double > row_max;
	for (int y = 0; y < h; y++)
		row_max.push_back(MaxElement(image[y]));

	return MaxElement(row_max);
}

double getMinValue( Buffer2d& image )
{
	int h = image.size();

	std::vector< double > row_min;
	for (int y = 0; y < h; y++)
		row_min.push_back(MinElement(image[y]));

	return MinElement(row_min);
}

// Regions
double maxValueInRegion( Buffer2d& image,  std::vector< Vec2i >& region )
{
	double result = 0;
	for (int i=0;i<region.size();i++)
	{
		int x = region[i].x();
		int y = region[i].y();

		if (image[y][x] > result)
			result = image[y][x];
	}

	return result;
}

Vec2i sizeofRegion( std::vector< Vec2i >& region )
{
	uint minX, maxX, minY, maxY;
	minX = maxX = region[0].x();
	minY = maxY = region[0].y();
	for (int i=1; i<region.size(); i++)
	{
		uint x = region[i].x();
		uint y = region[i].y();

		if (x < minX) minX = x;
		if (x > maxX) maxX = x;
		if (y < minY) minY = y;
		if (y > maxY) maxY = y;
	}

	return Vec2i( maxX - minX + 1, maxY - minY + 1 );
}

Vec2i centerOfRegion( std::vector< Vec2i >& region )
{
	Vec2i center(0, 0);

	if (!region.empty())
	{
		for (int i=0;i<region.size();i++)
			center += region[i];

		center /= region.size();
	}

	return center;
}

void BBofRegion( std::vector< Vec2i >& region, Vec2i &bbmin, Vec2i &bbmax )
{
	uint minX, maxX, minY, maxY;
	minX = maxX = region[0].x();
	minY = maxY = region[0].y();
	for (int i=1; i<region.size(); i++)
	{
		uint x = region[i].x();
		uint y = region[i].y();

		if (x < minX) minX = x;
		if (x > maxX) maxX = x;
		if (y < minY) minY = y;
		if (y > maxY) maxY = y;
	}

	bbmin = Vec2i(minX, minY);
	bbmax = Vec2i(maxX, maxY);
}

std::vector< double > getValuesInRegion( Buffer2d& image, 
	std::vector< Vec2i >& region, bool xFlipped /*= false*/ )
{
	std::vector< double > values;

	int w = image[0].size();
	uint x, y;
	for (int i = 0; i < region.size(); i++)
	{
		x = region[i].x();
		y = region[i].y();
		if(xFlipped) x = (w-1) - x;

		values.push_back(image[y][x]);
	}

	return values;
}

std::vector< Vec2i > getRegionGreaterThan( Buffer2d& image, Buffer2b& mask, Vec2i seed, double threshold )
{
	std::vector< Vec2i > region;

	int w = image[0].size();
	int h = image.size();		

	// Push the \seed to stack
	std::stack<Vec2i> activePnts;
	activePnts.push(seed);
	mask[seed.y()][seed.x()] = true;

	while (!activePnts.empty())
	{
		// Add one active point to \region
		Vec2i currP = activePnts.top();
		region.push_back(currP);
		activePnts.pop();

		// Push all its unvisited neighbors to the stack
		int min_x = RANGED(0, currP.x()-1 ,w-1);
		int max_x = RANGED(0, currP.x()+1, w-1);
		int min_y = RANGED(0, currP.y()-1 ,h-1);
		int max_y = RANGED(0, currP.y()+1, h-1);

		for (int y = min_y; y <= max_y; y++)
			for (int x = min_x; x <= max_x; x++)
			{
				if (!mask[y][x] && image[y][x]>threshold)
				{
					activePnts.push( Vec2i(x, y) );
					mask[y][x] = true;
				}
			}
	}

	return region;
}

std::vector< std::vector< Vec2i > > getRegionsGreaterThan( Buffer2d& image, double threshold )
{
	std::vector< std::vector< Vec2i > > regions;

	int w = image[0].size();
	int h = image.size();

	std::vector< std::vector< bool > > mask = createImage(w, h, false);

	for(int y = 0; y < h; y++){
		for(int x = 0; x < w; x++)	{
			if (!mask[y][x] && image[y][x]>threshold)
			{
				//saveAsImage(mask, "mask1.png");
				std::vector< Vec2i > region = getRegionGreaterThan(image, mask, Vec2i(x, y), threshold);
				regions.push_back(sampleRegion(region, 100));
				//regions.push_back(region);
				//saveAsImage(mask, "mask2.png");
			}
			mask[y][x] = true;
		}
	}

	// If there are a lot of hot regions, regard them as one
	std::vector< Vec2i > super_region, sampleSuper;
	if (regions.size() > 10)
	{
		foreach(std::vector< Vec2i > r, regions)
		{
			foreach(Vec2i p, r)
				super_region.push_back(p);
		}

		regions.clear();
		regions.push_back(sampleRegion(super_region, 100));
		//regions.push_back(super_region);
	}

	return regions;
}



// Shifting
std::vector< Vec2i > deltaVectorsToKRing( int deltaX, int deltaY, int K )
{
	std::vector< Vec2i > Vecs;

	int leftX = - deltaX * K;
	int rightX = -leftX;
	int topY = - deltaY * K;
	int bottomY = -topY;

	// top and bottom lines
	for (int x = leftX; x <= rightX; x += deltaX)
	{
		Vecs.push_back( Vec2i(x, topY) );
		Vecs.push_back( Vec2i(x, bottomY) );
	}

	// left and right lines
	for (int y = topY + deltaY; y <= bottomY - deltaY; y+= deltaY)
	{
		Vecs.push_back( Vec2i(leftX, y) );
		Vecs.push_back( Vec2i(rightX, y) );
	}

	return Vecs;
}

std::vector< Vec2i > shiftRegionInBB( std::vector< Vec2i >& region, Vec2i delta, Vec2i bbmin, Vec2i bbmax )
{
	std::vector< Vec2i > toRegion;

	for (int i = 0; i < region.size(); i++)
	{
		Vec2i p = region[i] + delta;
		if ( RANGE( p.x(), bbmin.x(), bbmax.x() ) && RANGE( p.y(), bbmin.y(), bbmax.y() ) )
			toRegion.push_back( p );
		else
		{// Out of BB
			toRegion.clear();
			break;
		}			
	}

	return toRegion;
}


// Color
void setRegionColor( Buffer2d& image, std::vector< Vec2i >& region, double color )
{
	uint w = image[0].size();
	uint h = image.size();

	for (int i = 0; i < region.size(); i++)
	{
		uint x = RANGED(0, region[i].x(), w-1);
		uint y = RANGED(0, region[i].y(), h-1);

		image[y][x] = color;
	}
}

void setPixelColor( Buffer2d& image, Vec2i pos, double color )
{
	uint w = image[0].size();
	uint h = image.size();

	uint x = RANGED(0, pos.x(), w-1);
	uint y = RANGED(0, pos.y(), h-1);

	image[y][x] = color;
}

QRgb jetColor( double val, double min, double max )
{
	uchar rgb[3];	

	ColorMap::jetColorMap(rgb, val, min, max);

	return QColor::fromRgb(rgb[0],rgb[1],rgb[2]).rgba();
}


// Visualize
void visualizeRegions( int w, int h, std::vector< std::vector<Vec2i> >& regions, QString filename )
{
	std::vector< std::vector< double > > debugImg = createImage(w, h, 0.0);
	double step = 1.0 / regions.size();
	for (int i=0;i<regions.size();i++)
	{
		setRegionColor(debugImg, regions[i], step * (i+1));
	}
	saveAsImage(debugImg, filename);
}

template< typename T >
std::vector< std::vector < T > > createImage( int w, int h, T intial )
{
	return std::vector< std::vector < T > > ( h, std::vector<T>( w, intial ) );
}

void saveAsImage( Buffer2d& image, QString fileName )
{
	int h = image.size();
	int w = image[0].size();
	QImage Output(w, h, QImage::Format_ARGB32);

	double minV = getMinValue(image);
	double maxV = getMaxValue(image);
	double range = maxV - minV;

	// \p y is flipped, since OpenGL has origin at the left bottom conner, while Qt at left top conner
	for(int y = 0; y < h; y++){
		for(int x = 0; x < w; x++)	{
			double c = (image[y][x] - minV) / range;
			Output.setPixel(x, (h-1)-y, jetColor( c, 0., 1.));			
		}
	}

	Output.save(fileName);
}

void saveAsData( Buffer2d& image, double maxV, QString fileName )
{
	int h = image.size();
	int w = image[0].size();

	QFile file(fileName); 
	file.open(QIODevice::WriteOnly | QIODevice::Text);

	// \p y is flipped, since OpenGL has origin at the left bottom conner, while Qt at left top conner
	for(int y = 0; y < h; y++){
		for(int x = 0; x < w; x++)	{
			file.write(qPrintable(QString::number(image[y][x]) + "\t"));
		}
		file.write("\n");
	}

	file.close();
}

void saveAsBinaryImage( Buffer2d& image, QString fileName )
{
	int h = image.size();
	int w = image[0].size();
	QImage Output(w, h, QImage::Format_ARGB32);

	QRgb red = QColor::fromRgb(255, 0, 0).rgba();
	QRgb blue = QColor::fromRgb(0, 0, 255).rgba();

	for(int y = 0; y < h; y++){
		for(int x = 0; x < w; x++)	{
			QRgb color = image[y][x] ? red : blue;
			Output.setPixel(x, (h-1)-y, color);			
		}
	}

	Output.save(fileName);
}

std::vector<int> uniformIntegerSamples( int N, int range )
{
	N = Min(N, range);

	std::vector<int> results;

	for (int i = 0; i < range; i++)
		results.push_back(i);

	// Shuffle
	std::random_shuffle ( results.begin(), results.end() );

	return std::vector<int>(results.begin(), results.begin() + N);
}

std::vector<Vec2i> sampleRegion( std::vector<Vec2i> &region, int N )
{
	if(region.size() <= N) return region;

	std::vector<Vec2i> samples;
	std::vector<int> sampleIdx = uniformIntegerSamples(N, region.size());
	foreach(int i, sampleIdx)
		samples.push_back(region[i]);

	return samples;
}



// Theta is measured in degree
Eigen::Matrix3d rotationMatrixAroundAxis( Vec3d u, double theta )
{
	u.normalize();

	double x = u[0], y = u[1], z = u[2];

	Eigen::Matrix3d I, cpm, tp, R;

	I = Eigen::Matrix3d::Identity(3,3);

	tp <<	x*x, x*y, x*z,
		x*y, y*y, y*z,
		x*z, y*z, z*z;

	cpm <<   0, -z,  y,
		z,  0, -x,
		-y,  x,  0;

	theta = RADIANS(theta);
	R = cos(theta)*I + sin(theta)*cpm + (1-cos(theta))*tp;

	return R;
}

Vec3d rotatePointByMatrix( Eigen::Matrix3d &R, Vec3d p )
{
	Eigen::Vector3d rp = R * V2E(p);
	return E2V(rp);
}

double volumeOfBB( Vec3d &extents )
{
	return extents.x() * extents.y() * extents.z();
}


void twoFurthestPoints( std::vector<Point> &points, Point &p1, Point &p2 )
{
	if(points.size() < 2) return;

	double maxDis = 0.0;
	int id1, id2;
	for (int i = 0; i < points.size(); i++)
	{
		for (int j = i + 1; j <points.size(); j++)
		{
			double dis = (points[i]-points[j]).norm();
			if (dis > maxDis)
			{
				id1 = i;
				id2 = j;
				maxDis = dis;
			}
		}
	}

	p1 = points[id1];
	p2 = points[id2];
}

double distanceCluster2Cluster( std::vector<Point> &cluster1, std::vector<Point> &cluster2 )
{
	double minDis = DOUBLE_INFINITY;
	foreach(Point p1, cluster1)
	{
		foreach(Point p2, cluster2)
		{
			double dis = (p1-p2).norm();

			if (dis < minDis)
			{
				minDis = dis;
			}
		}
	}

	return minDis;
}

std::vector< std::vector<Point> > twoClassClustering( std::vector<Point>& points, Point seed1, Point seed2 )
{
	// Assign all the points to two seeds
	std::vector<  std::vector< Point > >clusters(2);
	foreach( Vec3d p, points)
	{
		if ((p - seed1).norm() < (p - seed2).norm())
			clusters[0].push_back(p);
		else
			clusters[1].push_back(p);
	}

	return clusters;
}

Point centroid(std::vector<Point> points)
{
	Point center(0);
	foreach(Point p, points) center += p;
	center /= points.size();

	return center;
}

std::vector<Point> uniformSampleCurve(std::vector<Point> & points)
{
	std::vector<Point> newPoints;

	// Cumulative length
	int N = points.size();
	std::vector<double> cumDis(N, 0.0);
	for (int i = 1; i < N; i++)
	{
		double len = (points[i] - points[i-1]).norm();
		cumDis[i] = cumDis[i-1] + len;
	}

	// Sampling
	int low = 0; // the starting end
	double step = cumDis[N-1] / (N-1);
	newPoints.push_back(points[0]);
	double pos = 0; // current length
	for (int i = 1; i < N; i++)
	{
		pos += step;
		while(low < N-2 && cumDis[low+1] < pos)
			low ++;

		double a = pos - cumDis[low];
		double b = cumDis[low+1] - pos;
		double alpha = a / (a+b);
		double belta = b / (a+b);

		Point p = belta * points[low] + alpha * points[low + 1];
		newPoints.push_back(p);
	}

	return newPoints;
}

Vec3d maxMidMinValues( Buffer2d & image )
{
	//std::set<double> vals;

	//int h = image.size();
	//int w = image[0].size();

	//for(int y = 0; y < h; y++){
	//	for(int x = 0; x < w; x++)	{
	//		vals.insert(image);
	//	}
	//}

	//if(vals.size() > 1)
	//	return Vec3d(*vals.begin(), *(vals.begin()++), *vals.rbegin());
	//else
	//	return Vec3d(*vals.begin(),0,0);
	return Vec3d(0);
}

Buffer2v2i getMaximumRegions( Buffer2d &image )
{
	// Precondition: each pixel is greater or equal than 0
	double maxV = getMaxValue(image);
	Buffer2v2i regions;
	double hot_cap = 1.0;

	while (regions.empty())
	{
		hot_cap -= 0.05; // increase the cap
		regions = getRegionsGreaterThan(image, maxV * hot_cap);

		// If the hot regions are too small
		int num = 0;
		for (int i = 0; i < regions.size(); i++)
		{
			if (regions[i].size() > num) num = regions[i].size();
		}

		if(num < 3) 
			regions.clear();	
	}

	return regions;
}


// AABB
void computeAABB( std::vector<Point> &points,  Vec3d &bbmin, Vec3d &bbmax)
{
	bbmin = Point( FLT_MAX,  FLT_MAX,  FLT_MAX);
	bbmax = Point(-FLT_MAX, -FLT_MAX, -FLT_MAX);	

	foreach(Point p, points)
	{
		bbmin.minimize(p);
		bbmax.maximize(p);
	}
}

std::vector<Point> cornersOfAABB(Vec3d bbmin, Vec3d bbmax)
{
	std::vector<Point> corner;

	double min_x = bbmin.x();
	double min_y = bbmin.y();
	double min_z = bbmin.z();	
	
	double max_x = bbmax.x();
	double max_y = bbmax.y();
	double max_z = bbmax.z();

	corner.push_back(Point( min_x, min_y, min_z));
	corner.push_back(Point( min_x, min_y, max_z));
	corner.push_back(Point( min_x, max_y, min_z));
	corner.push_back(Point( min_x, max_y, max_z));

	corner.push_back(Point( max_x, min_y, min_z));
	corner.push_back(Point( max_x, min_y, max_z));
	corner.push_back(Point( max_x, max_y, min_z));
	corner.push_back(Point( max_x, max_y, max_z));

	return corner;
}