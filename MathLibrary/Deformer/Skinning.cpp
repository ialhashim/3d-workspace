#include "Skinning.h"
#include "GraphicsLibrary/Basic/Line.h"
#include "DualQuat.h"
#include "Utility/Macros.h"
#include "GraphicsLibrary/Basic/Plane.h"

Skinning::Skinning( QSurfaceMesh * src_mesh, GeneralizedCylinder * using_gc )
{
	this->mesh = src_mesh;
	this->currGC = using_gc;
	this->origGC = *currGC;

	computeMeshCoordinates();
}

Skinning::SkinningCoord Skinning::computeCoordinates( GeneralizedCylinder *gc, Point& v )
{
	// Go over the skeleton for the closest segment
	double minDis = DBL_MAX, minT = 0;
	int minIdx = -1;
	for(int i = 0; i < (int)gc->crossSection.size() - 1; i++)
	{
		GeneralizedCylinder::Circle & c1 = gc->crossSection[i];
		GeneralizedCylinder::Circle & c2 = gc->crossSection[i+1];

		Line seg(c1.center, c2.center);
		double t = seg.timeAt(v);

		// Check belongs to segment
		if(t >= 0.0 && t <= 1.0){

			double dist = seg.distanceToUnbounded(v);
			double radAt = ((c1.radius * (1-t)) + (c2.radius * t)) * 1.2;

			// Check inside GC
			if( dist < radAt && dist < minDis )
			{
				minDis = dist;
				minIdx = i;
				minT = t;
			}
		}
	}

	if(minIdx == -1)
	{
		// At the end of or outside \gc
		// Find the closest cross section
		minDis = DBL_MAX;
		for(int j = 0; j < (int)gc->crossSection.size(); j++)
		{
			double dist = (gc->crossSection[j].center - v).norm();

			if(dist < minDis){
				minDis = dist;
				minIdx = j;
			}
		}

		GeneralizedCylinder::Circle & c = gc->crossSection[minIdx];
		Line seg(c.center, c.center + c.normal());
		double t = seg.timeAt(v);
		Point proj = c.center + c.normal()*t;

		return SkinningCoord(minIdx, minIdx, t, v - proj);
	}
	else
	{
		// On the segment [\minIdx, \minIdx+1]
		Point proj = gc->crossSection[minIdx].center*(1-minT) + gc->crossSection[minIdx+1].center*minT;
		return SkinningCoord(minIdx, minIdx + 1, minT, v - proj);
	}
}

void Skinning::computeMeshCoordinates()
{
	// The coordinates are computed based the original GC \origGC
	Surface_mesh::Vertex_property<Point> points = mesh->vertex_property<Point>("v:point");
	Surface_mesh::Vertex_iterator vit, vend = mesh->vertices_end();

	coordinates.clear();
	for (vit = mesh->vertices_begin(); vit != vend; ++vit)
	{
		Vec3d v = points[vit];
		coordinates.push_back(computeCoordinates(&origGC, v));
	}
}

Point Skinning::fromCoordinates( GeneralizedCylinder &orig_gc, SkinningCoord coords )
{
	int i1 = coords.n1;
	int i2 = coords.n2;
	double w = coords.time;

	GeneralizedCylinder::Circle & orig_c1 = orig_gc.crossSection[i1];
	GeneralizedCylinder::Circle & orig_c2 = orig_gc.crossSection[i2];
	GeneralizedCylinder::Circle & c1 = currGC->crossSection[i1];
	GeneralizedCylinder::Circle & c2 = currGC->crossSection[i2];

	// Projection on the skeleton
	Point proj;
	if (i1 == i2)
		proj = orig_c1.center + orig_c1.normal() * w;
	else
	{
		Line seg(orig_c1.center, orig_c2.center);
		proj = seg.pointAt(w);
	}

	// Scaling
	double s;
	if(i1 == i2)
		s = c1.radius / orig_c1.radius;
	else
	{
		// Interpolate radius
		double orig_r = orig_c1.radius*(1.0-w) + orig_c2.radius*w;
		double r = c1.radius*(1.0-w) + c2.radius*w;	 
		s = r / orig_r;
	}

	// Rotation and translation 
	// Using dual quaternion blending
	Matrix3d R1 = rotationOfCurve(i1);
	Matrix3d R2 = rotationOfCurve(i2);
	Vector3d T1 = V2E(c1.center) - R1 * V2E(orig_c1.center);
	Vector3d T2 = V2E(c2.center) - R2 * V2E(orig_c2.center);

	DualQuat dq1, dq2, dqb;
	dq1.SetTransform(R1, T1);
	dq2.SetTransform(R2, T2);
	dqb = dq1*(1-w)  +  dq2*w;

	Matrix3d R; Vector3d T;
	dqb.GetTransform(R, T); 

	// The result
	Vector3d V = V2E((proj + coords.d * s));
	V = R * V + T;

	Vec3d point = E2V(V);
	if (!( point[0] == point[0]))
	{
		int aaaa = 0;
	}


	return point;
}

std::vector<double> Skinning::getCoordinate( Point p )
{
	// The coordinates in the \currGC
	std::vector<double> coords;

	SkinningCoord skinning_coord = computeCoordinates(currGC, p);

	coords.push_back(skinning_coord.n1);
	coords.push_back(skinning_coord.n2);
	coords.push_back(skinning_coord.time);

	coords.push_back(skinning_coord.d[0]);
	coords.push_back(skinning_coord.d[1]);
	coords.push_back(skinning_coord.d[2]);

	// The \currGC has also to be stored as part of the coordinates
	// It will be used to compute the relative local rotation and translation
	foreach (GeneralizedCylinder::Circle c, currGC->crossSection )
	{
		// center
		coords.push_back(c.center[0]);
		coords.push_back(c.center[1]);
		coords.push_back(c.center[2]);

		// normal
		coords.push_back(c.n[0]);
		coords.push_back(c.n[1]);
		coords.push_back(c.n[2]);

		// radius
		coords.push_back(c.radius);
	}

	return coords;
}

Point Skinning::fromCoordinates( std::vector<double> &coords )
{
	int N = currGC->crossSection.size();

	// 6 for \SkinningCoord
	// 7 for each \Circle of cross sections
	if(coords.size() != 6+7*N ) return Point(0.0);

	// The \SkinningCoord
	SkinningCoord skinning_coords;
	skinning_coords.n1 = int(coords[0]);
	skinning_coords.n2 = int(coords[1]);
	skinning_coords.time = coords[2];
	skinning_coords.d = Vec3d(coords[3], coords[4], coords[5]);

	// The original GC
	GeneralizedCylinder orig_gc = *currGC;
	int id = 6;
	for (int i = 0; i < N; i++, id += 7)
	{
		GeneralizedCylinder::Circle &c = orig_gc.crossSection[i];
		c.center = Point(coords[id], coords[id+1], coords[id+2]);
		c.n = Vec3d(coords[id+3], coords[id+4], coords[id+5]);
		c.radius = coords[id+6];
	}

	return fromCoordinates(orig_gc, skinning_coords);
}

Matrix3d Skinning::rotationOfCurve( int cid )
{
	Matrix3d rm = Matrix3d::Identity();

	if (cid < 0 || cid >= currGC->crossSection.size())
		return rm;

	GeneralizedCylinder::Circle & orig_c = origGC.crossSection[cid];
	GeneralizedCylinder::Circle & c = currGC->crossSection[cid];

	Vec3d n1(orig_c.normal()), n2(c.normal());
	Vector3d axis = V2E(cross(n1, n2).normalized());

	if(axis.norm() > 0)
	{
		double angle = acos(RANGED(-1.0, dot(n1, n2), 1.0));
		rm = AngleAxisd(angle, axis).toRotationMatrix();
	}

	return rm;
}

void Skinning::deform()
{
	Surface_mesh::Vertex_property<Point> points = mesh->vertex_property<Point>("v:point");
	Surface_mesh::Vertex_iterator vit, vend = mesh->vertices_end();

	for (vit = mesh->vertices_begin(); vit != vend; ++vit)
	{		
		int vi = Surface_mesh::Vertex(vit).idx();
		points[vit] = fromCoordinates(origGC, coordinates[vi]);
	}
}

bool Skinning::atEnd( Point p )
{
	SkinningCoord skinning_coord = computeCoordinates(currGC, p);
	int n = skinning_coord.n1;
	int N = currGC->frames.count() - 1;
	double pos = (double)n / (double)N;
	if (pos < 0.2 || pos > 0.8 || n < 2 || n > N-1)
		return true;
	else
		return false;
}

Point Skinning::closestProjection( Point v )
{
	SkinningCoord coords = computeCoordinates(currGC, v);

	int i1 = coords.n1;
	int i2 = coords.n2;
	double w = coords.time;
	Vec3d d = coords.d;

	GeneralizedCylinder::Circle & c1 = currGC->crossSection[i1];
	GeneralizedCylinder::Circle & c2 = currGC->crossSection[i2];

	// Blended Radius
	double r = c1.radius*(1.0-w) + c2.radius*w;

	// The Projection
	Point p;
	if ( (i1 == 0 && i2 == 0)
		|| i1 == currGC->crossSection.size()-1 )
	{
		// Out the ends of gc
		// Project within the circle
		Vec3d n = currGC->crossSection[i1].n;
		Point c = currGC->crossSection[i1].center;
		Plane plane(n, c);
		Point proj = plane.projectionOf(v);
		Vec3d delta = c - proj;
		p = Min(r, delta.norm()) * delta.normalized();
	}
	else
	{
		// Within two ends
		// Projection on the skeleton
		Point skl_proj;
		if (i1 == i2)
			skl_proj = c1.center + c1.normal() * w;
		else
		{
			Line seg(c1.center, c2.center);
			skl_proj = seg.pointAt(w);
		}

		// Project to the GC surface
		p = skl_proj + r * d.normalized();
	}

	return p;
}
