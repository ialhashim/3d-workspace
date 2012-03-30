#include "Triangle.h"

float Triangle::edgeLenY(int e) const
{
	if(e == 0)			return p[1].y() - p[0].y();
	else if(e == 1)		return p[2].y() - p[1].y();
	else				return p[0].y() - p[2].y();
}

void Triangle::unloop(float threshold, Triangle & other)
{
	int edge1 = edgeLenY(0);
	int edge2 = edgeLenY(1);
	int edge3 = edgeLenY(2);

	// Trivial case: non looping triangle
	if(edge1 < threshold && edge2 < threshold && edge3 < threshold)
		return;

	float fullHeight = threshold * 2.0f;

	if(abs(edge1) < threshold)
	{
		// then 'p3' is looped
		if(edge3 > 0)	p[2].y() += fullHeight;
		else p[2].y() -= fullHeight;

		if(p[2].y() > 0)	other = this->shiftY(-fullHeight);
		else			other = this->shiftY(fullHeight);
	}

	if(abs(edge2) < threshold)
	{
		// then 'p1' is looped
		if(edge1 > 0)	p[0].y() += fullHeight;
		else p[0].y() -= fullHeight;

		if(p[0].y() > 0)	other = this->shiftY(-fullHeight);
		else			other = this->shiftY(fullHeight);
	}

	if(abs(edge3) < threshold)
	{
		// then 'p2' is looped
		if(edge2 > 0)	p[1].y() += fullHeight;
		else p[1].y() -= fullHeight;

		if(p[1].y() > 0)	other = this->shiftY(-fullHeight);
		else			other = this->shiftY(fullHeight);
	}

	other.flag = MODIFIED_TRI;
}

Triangle Triangle::shiftY( float offsetY )
{
	Triangle temp = *this;

	temp.p[0].y() += offsetY;
	temp.p[1].y() += offsetY;
	temp.p[2].y() += offsetY;

	return temp;
}

void Triangle::intersectionTest(const Ray & ray, HitResult & res, bool allowBack) const
{
	res.hit = false;
    res.distance = DBL_MAX;

	Vec3d vertex1 = p[0];
	Vec3d vertex2 = p[1];
	Vec3d vertex3 = p[2];

	// Compute vectors along two edges of the triangle.
	Vec3d edge1 = vertex2 - vertex1;
	Vec3d	edge2 = vertex3 - vertex1;

	// Compute the determinant.
	Vec3d directionCrossEdge2 = cross(ray.direction, edge2);

	float determinant = dot(edge1, directionCrossEdge2);

	// If the ray is parallel to the triangle plane, there is no collision.
	if (abs(determinant) < 1e-7)
		return;

	float inverseDeterminant = 1.0f / determinant;

	// Calculate the U parameter of the intersection point.
	Vec3d distanceVec3dtor = ray.origin - vertex1;
	float triangleU = dot(distanceVec3dtor, directionCrossEdge2);
	triangleU *= inverseDeterminant;

	// Make sure it is inside the triangle.
	if (triangleU < 0 || triangleU > 1)
		return;

	// Calculate the V parameter of the intersection point.
	Vec3d distanceCrossEdge1 = cross(distanceVec3dtor, edge1);
	float triangleV = dot(ray.direction, distanceCrossEdge1);
	triangleV *= inverseDeterminant;

	// Make sure it is inside the triangle.
	if (triangleV < 0 || triangleU + triangleV > 1)
		return;

	// Compute the distance along the ray to the triangle.
	float rayDistance = dot(edge2, distanceCrossEdge1);
	rayDistance *= inverseDeterminant;

	if(!allowBack)
	{
		// Is the triangle behind the ray origin?
		if (rayDistance < 0)
			return;
	}

	res.hit = true;
	res.distance = rayDistance;

	res.u = triangleU;
	res.v = triangleV;

	res.index = this->index;
}
