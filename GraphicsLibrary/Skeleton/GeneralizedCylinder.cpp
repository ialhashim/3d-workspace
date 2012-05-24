#include "GeneralizedCylinder.h"
#include "RMF.h"
#include "GraphicsLibrary/SpacePartition/Intersection.h"
#include "ClosedPolygon.h"
#include "Utility/SimpleDraw.h"

GeneralizedCylinder::GeneralizedCylinder( std::vector<Point> spinePoints, QSurfaceMesh * mesh, bool computeRadius /*= true */ )
{
	// Build minimum rotation frames on this spine
	frames = RMF (spinePoints);

	double sumRadius = 0; // for average radius
	int numNonZero = 0;
	double nonZeroRadius = 1.0;

	// Build cross-section
	for(uint i = 0; i < spinePoints.size(); i++)
	{
		// Find radius from the cross-section
		double radius = 0;

		if(computeRadius)
		{
			std::vector<Point> cur_cs;
		
			ClosedPolygon polygon(frames.point[i]);

			for(uint fi = 0; fi < mesh->face_array.size(); fi++)
			{
				Vec3d p1, p2;

				if(ContourFacet(frames.U[i].t, frames.point[i], mesh->facePoints(mesh->face_array[fi]), p1, p2) > 0)
					polygon.insertLine(p1,p2);
			}

			polygon.close();

			// Sort based on distance and filter based on segment distance
			foreach(Point p, polygon.closedPoints)
			{
				radius = Max(radius, (p - frames.point[i]).norm());
			}

			// Some filters to be fail-safe
			if(nonZeroRadius == 0 && radius != 0)	nonZeroRadius = radius;
			if(radius > 10 * nonZeroRadius)	radius = 0;
		}

		crossSection.push_back(Circle(frames.point[i], radius, frames.U[i].t, i));

		sumRadius += radius;
		if(radius != 0) numNonZero++;
	}

	// Start and end
	if(crossSection.front().radius == 0) 
		crossSection.front().radius = crossSection[1].radius;

	if(crossSection.back().radius == 0) 
		crossSection.back().radius = crossSection[crossSection.size() - 2].radius;

	// Some filters to be fail-safe
	double avgRadius = sumRadius / numNonZero;
	for (uint i = 0; i < crossSection.size(); i++)
		if(crossSection[i].radius == 0) crossSection[i].radius = avgRadius;

	// Smoothing for some robustness
	std::vector<double> filterRadius(crossSection.size());
	for(int k = 0; k < 2; k++)
	{
		for (uint i = 0; i < crossSection.size(); i++)
			filterRadius[i] = crossSection[i].radius;

		for (uint i = 1; i < crossSection.size() - 1; i++)
			filterRadius[i] = (filterRadius[i-1] + filterRadius[i+1]) / 2.0;
		filterRadius.front() = filterRadius[0]; filterRadius.back() = filterRadius[filterRadius.size() - 2];

		for (uint i = 0; i < crossSection.size(); i++)
			crossSection[i].radius = filterRadius[i];
	}

	// push ends a tiny bit
	crossSection.front().center -= 0.01 * (crossSection[1].center - crossSection[0].center).normalized();
	crossSection.back().center += 0.01 * (crossSection.back().center - crossSection[crossSection.size() - 2].center).normalized();

	src_mesh = mesh;
	isDrawFrames = false;

	printf("\nGC generated.\n");
}

void GeneralizedCylinder::draw()
{
	if(!frames.U.size()) return;

	glDisable(GL_DEPTH_TEST);

	foreach(Circle c, crossSection)
	{
		std::vector<Point> pnts = c.toSegments(30, frames.U[c.index].s);
		pnts.push_back(pnts.front());

		SimpleDraw::IdentifyConnectedPoints(pnts, 0,0,0);
	}

	foreach(Point p, debugPoints)
		SimpleDraw::IdentifyPoint(p);

	// Draw frames
	if(isDrawFrames)
	{
		std::vector<Point> dir1, dir2, dir3;
		for(uint i = 0; i < frames.U.size(); i++){
			dir1.push_back(frames.U[i].r);
			dir2.push_back(frames.U[i].s);
			dir3.push_back(frames.U[i].t);
		}
		SimpleDraw::DrawLineTick(frames.point, dir1, 0.1f, false, 1,0,0,1);
		SimpleDraw::DrawLineTick(frames.point, dir2, 0.1f, false, 0,1,0,1);
		SimpleDraw::DrawLineTick(frames.point, dir3, 0.1f, false, 0,0,1,1);
	}

	glEnable(GL_DEPTH_TEST);
}

std::vector<Point> GeneralizedCylinder::Circle::toSegments( int numSegments, const Vec3d& startVec, double delta )
{
	std::vector<Point> result;
	double theta = (2 * M_PI) / numSegments;
	Vec3d vi = startVec;

	for(int i = 0; i < numSegments; i++)
	{
		result.push_back(center + (vi * radius * delta) );

		// Rodrigues' rotation formula
		vi = vi * cos(theta) + cross(n, vi) * sin(theta) + n * dot(n, vi) * (1 - cos(theta));
	}

	return result;
}

Normal GeneralizedCylinder::Circle::normal()
{
	return n;
}

void GeneralizedCylinder::realignCrossSections()
{
	for(uint i = 0; i < crossSection.size(); i++)
	{
		crossSection[i].center = frames.point[i];
		crossSection[i].n = frames.U[i].t;
	}
}

