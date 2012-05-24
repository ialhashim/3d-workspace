#include "ClosedPolygon.h"
#include "GraphicsLibrary/Basic/PolygonArea.h"
#include "Utility/Macros.h"


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

void ClosedPolygon::computeLengths()
{
	int N = closedPoints.size();

	closedLength = 0;

	for(int i = 0; i < N; i++)
	{
		double dist = (closedPoints[(i+1) % N] - closedPoints[i]).norm();

		minEdgeLength = Min(minEdgeLength, dist);

		closedEdgeLen.push_back(dist);
		closedLength += dist;
	}
}

std::vector<Vec3d> ClosedPolygon::getEqualDistancePoints(int numSides, const Vec3d& center)
{
	std::vector<Vec3d> result;

	int N = closedPoints.size();

	if(N < 1)	return result; // empty polygon

	for(int i = 0; i < N; i++)
		lines.push_back(Line(closedPoints[i], closedPoints[(i+1) % N], i));

	this->computeLengths();

	// Distance to walk on polygon
	double segmentLength = this->closedLength / numSides;

	// Locate start point using vecUp
	Vec3d startPoint;
	int startIndex = 0;

	Plane halfPlane(vecUp, center);
	//testPlanes1.push_back(halfPlane);

	double minDist = DBL_MAX;

	// Test intersection with all lines and remember minimum one
	for(int i = 0; i < N; i++)
	{
		Vec3d pointIntersect;

		int res = halfPlane.LineIntersect(lines[i], pointIntersect);

		if(res == INTERSECT || res == ENDPOINT_INTERSECT)
		{
			Vec3d toIntsect = pointIntersect - center;

			if(toIntsect.norm() < minDist && dot(toIntsect, vecB) > 0)
			{
				minDist = toIntsect.norm();

				startPoint = pointIntersect;
				startIndex = i;
			}
		}
	}

	double t = lines[startIndex].timeAt(startPoint);
	int index = startIndex;

	// Compute equal-dist points on polygon
	for(int s = 0; s < numSides; s++)
	{
		// Add new point
		result.push_back(lines[index].pointAt(t));

		walk(segmentLength, t, index, &t, &index);
	}

	// if polygon is opposite direction then reverse 
	if( signedArea(result, plane.n, center) < 0 )
	{
		std::reverse(result.begin(), result.end());
		std::rotate(result.begin(), result.begin()+result.size()-1 , result.end());
	}

	return closedPoints = result;
}

void ClosedPolygon::walk(double distance, double startTime, int index, double * destTime, int * destIndex)
{
	double remain = lines[index].lengthsAt(startTime).second;

	// Case 1: the point is on the starting line
	if(remain > distance)
	{
		double startLength = startTime * lines[index].length;
		*destTime = (startLength + distance) / lines[index].length;
		*destIndex = index;
		return;
	}

	double walked = remain;

	// Case 2: keep walking next lines
	while(walked < distance)
	{
		index = (index + 1) % lines.size();		// step to next line
		walked += lines[index].length;
	}

	// Step back to the start of this line
	walked -= lines[index].length;

	double remainDistance = distance - walked;
	double endTime = remainDistance / lines[index].length;

	*destTime = endTime;
	*destIndex = index;
}
