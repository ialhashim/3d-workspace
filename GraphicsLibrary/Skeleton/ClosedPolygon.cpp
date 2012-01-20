#include "ClosedPolygon.h"

ClosedPolygon::ClosedPolygon( const Vec3d& centerPoint, double epsilon)
{
	this->points.create(3);

	lastVertexIndex = lastEdgeIndex = 0;

	center = centerPoint;

	closedPolyEpsilon = epsilon;
}

void ClosedPolygon::insertPoint(const Vec3d & point)
{
	closedPoints.push_back(point);
}

int ClosedPolygon::insertIndexedPoint(const Vec3d &p)
{
	int vIndex = -1;

	kdres * findP = points.nearest3f(p.x(), p.y(), p.z());

	if(findP)
	{
		double * pos = findP->riter->item->pos;

		Vec3d closestPoint(pos[0], pos[1], pos[2]);

		double dist = (closestPoint - p).norm();

		if(dist < closedPolyEpsilon)
		{
			vIndex = findP->riter->item->index;
			kd_res_free(findP);
			return vIndex;
		}
	}

	vIndex = lastVertexIndex;
	allPoints[vIndex] = p;

	points.insert3f(p.x(), p.y(), p.z(), lastVertexIndex++);

	return vIndex;
}

void ClosedPolygon::insertLine(const Vec3d &p1, const Vec3d &p2)
{
	int vIndex1 = insertIndexedPoint(p1);
	int vIndex2 = insertIndexedPoint(p2);

	g.AddEdge(vIndex1, vIndex2, 1, lastEdgeIndex++);
}

bool ClosedPolygon::isClosed(int& i)
{
	//i = g.getNodeLargestConnected();

	double minDist = DBL_MAX;

	for(std::map<int, Vec3d>::iterator pnt = allPoints.begin(); pnt != allPoints.end(); pnt++){
		double curDist = (center - pnt->second).norm();

		if(curDist < minDist){
			minDist = curDist;
			i = pnt->first;
		}
	}

	if(i == -1) return false;

	return g.isCircular(i);
}

void ClosedPolygon::close()
{
	if(!allPoints.size())
		return;

	int i = 0;

	// for non circular cross-sections, force closed circle (won't work with multi-non connected)
	if(!isClosed(i))
	{
		std::vector<uint> leaves = g.GetLeaves();

		// connect ends
		if(leaves.size() > 1)
			g.AddEdge(leaves.front(), leaves.back(), 1, lastEdgeIndex++);
	}

	if(i >= 0)
	{
		int j = g.GetRandomNeighbour(i);

		g.SetEdgeWeight(i, j, DBL_MAX);

		std::list<uint> points = g.DijkstraShortestPath(i,j);

		for(std::list<uint>::iterator it = points.begin(); it != points.end(); it++)
			closedPoints.push_back(allPoints[*it]);
	}
}
