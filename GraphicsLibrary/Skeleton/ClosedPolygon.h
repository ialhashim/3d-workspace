#include "Vector.h"
#include "Utility/Graph.h"
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

	KDTree points;
	Vec3d center;

	Graph<uint, double> g;
	int lastVertexIndex;
	int lastEdgeIndex;

	std::map<int, Vec3d> allPoints;
	std::vector<Vec3d> closedPoints;
	std::vector<Vec3d> unvisitedPoints;

	std::vector<Vec3d> debugPoints;
};
