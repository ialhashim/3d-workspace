#include "Box3.h"

#include <algorithm>

Box3::Box3( Box3 &box )
{
	*this = box;
}

Box3::Box3()
{
	Axis.resize(3);
}

bool Box3::operator==( Box3& box )
{
	return (Center == box.Center)
		&& (Axis[0] == box.Axis[0])
		&& (Axis[1] == box.Axis[1])
		&& (Axis[2] == box.Axis[2])
		&& (Extent[0] == box.Extent[0])
		&& (Extent[1] == box.Extent[1])
		&& (Extent[2] == box.Extent[2]);
}

Vec3d Box3::ClosestPoint( const Vec3d& p )
{
	Vec3d d = p - Center;

	// Start result at center of box; make steps from there
	Vec3d q = Center;

	// For each OBB axis...
	for (int i = 0; i < 3; i++) {
		// ...project d onto that axis to get the distance
		// along the axis of d from the box center
		double dist = dot(d, Axis[i]);

		// If distance farther than the box extents, clamp to the box
		if (dist > Extent[i]) dist = Extent[i];
		if (dist < -Extent[i]) dist = -Extent[i];

		// Step that distance along the axis to get world coordinate
		q += dist * Axis[i];
	}

	return q;
}

void Box3::ClosestSegment( Box3 other, Vec3d & p, Vec3d & q)
{
	p = ClosestPoint(other.Center);
	q = other.ClosestPoint(p);
	p = ClosestPoint(q);
}



void Box3::normalizeAxis()
{
	for (int i=0;i<3;i++)	
		Axis[i].normalize();
}

void Box3::makeRightHanded()
{
	if ( dot( Axis[2], cross(Axis[0], Axis[1]) ) < 0 )
		Axis[2] *= -1;
}

// In ascending order
void Box3::sort()
{
	// Sort extents
	std::vector<double> ext;
	ext.push_back(Extent[0]);
	ext.push_back(Extent[1]);
	ext.push_back(Extent[2]);
	
	std::vector<double> sortedExt = ext;
	std::sort(sortedExt.begin(), sortedExt.end());

	Extent = Vec3d(sortedExt[0], sortedExt[1], sortedExt[2]);

	// Sorting indices
	std::vector<int> indices;
	foreach(double e, sortedExt)
		indices.push_back(std::find(ext.begin(), ext.end(), e) - ext.begin());

	// Sort axis
	std::vector<Vec3d> axisCopy = Axis;
	Axis.clear();
	foreach(int id, indices)
		Axis.push_back(axisCopy[id]);
}

Vec3d Box3::ClosestAxis( const Vec3d& v )
{
	double maxDot = 0;
	int bestID = 0;
	for (int i = 0; i< 3; i++)
	{
		double dp = abs(dot(Axis[i], v));
		if (dp > maxDot)
		{
			maxDot = dp;
			bestID = i;
		}
	}

	return Axis[bestID];
}
