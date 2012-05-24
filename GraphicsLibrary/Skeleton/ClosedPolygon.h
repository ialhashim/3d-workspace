#include "Vector.h"
#include "Utility/Graph.h"
#include "GraphicsLibrary/Basic/Line.h"
#include "GraphicsLibrary/Basic/Plane.h"
#include "GraphicsLibrary/SpacePartition/kdtree.h"
typedef unsigned int uint;

class ClosedPolygon{
public:
	ClosedPolygon( const Vec3d& centerPoint, double epsilon = 1.0e-7f);

	double closedPolyEpsilon;

	int insertIndexedPoint(const Vec3d &);
	void insertPoint(const Vec3d & point);
	void insertLine(const Vec3d &, const Vec3d &);
	
	void close();
	bool isClosed(int& i);

	std::vector<Vec3d> getEqualDistancePoints(int numSides, const Vec3d& center);
	void walk(double distance, double startTime, int index, double * destTime, int * destIndex);
	void computeLengths();
	std::vector<double> closedEdgeLen;
	double closedLength;
	double minEdgeLength;

	KDTree points;
	Vec3d center;

	Graph<uint, double> g;
	int lastVertexIndex;
	int lastEdgeIndex;

	std::map<int, Vec3d> allPoints;
	Plane plane;
	Vec3d vecUp, vecB;

	std::vector<Line> lines,tempLines,allLines;
	std::vector<Vec3d> closedPoints;
	std::vector<Vec3d> unvisitedPoints;

	std::vector<Vec3d> debugPoints;
};
