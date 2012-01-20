#include "VoxelDeformer.h"
#include "SimpleDraw.h"
#include "ColorMap.h"

VoxelDeformer::VoxelDeformer( QSurfaceMesh * fromMesh, double voxel_size )
{
	this->mesh = fromMesh;
	this->voxelSize = voxel_size;
	this->voxelRange = 0;

	this->voxler = new Voxeler(mesh, voxelSize, true);

	// Grow by specified range
	for(int i = 0; i < voxelRange; i++)
		this->voxler->grow();
	this->voxler->update();

	uint NV = voxler->voxels.size();

	meshPoints = mesh->clonePoints();
	pointToFFD.resize(meshPoints.size());

	QMap< int, std::map<int, Point> > pnts;

	// Find mesh points for each voxel
	for(int h = 0; h < (int) meshPoints.size(); h++ )
	{
		Point point = meshPoints[h];

		std::map<int, Voxel> neigh = voxler->around(point);
		for(std::map<int, Voxel>::iterator it = neigh.begin(); it!=neigh.end();it++)
		{
			int vi = it->first;
			
			Voxel v = voxler->voxels[vi];
			Point center(v.x * voxelSize, v.y * voxelSize, v.z * voxelSize);

			bool isInside = false;

			if(abs(center.x() - point.x()) < voxelSize * 0.5 && 
				abs(center.y() - point.y()) < voxelSize * 0.5 && 
				abs(center.z() - point.z()) < voxelSize * 0.5)
			{
				isInside = true;
				pnts[vi][h] = point;
				pointToFFD[h] = vi;
				break;
			}
		}
	}
	
	// build FFD cages
	for(int i = 0; i < (int)voxler->voxels.size(); i++)
	{
		Voxel v = voxler->voxels[i];

		Vec3d pos(v.x, v.y, v.z); pos *= voxelSize;

		ffd.push_back(new FFD());
		ffd.back()->fixed(Vec3i(2,2,2), pos, voxelSize, pnts[i]);
	}
	
	voxelPnts = voxler->corners;

	// Find correspondence
	for(uint i = 0; i < voxelPnts.size(); i++)
	{
		QControlPoint * closestPoint;
		Point corner = voxelPnts[i];
		cpnts.push_back(new QControlPoint(corner, i, Vec3i(0,0,0), 0));
		connect(cpnts.back(), SIGNAL(manipulated()), SLOT(update()));

		for(StdVector<FFD*>::iterator f = ffd.begin(); f != ffd.end(); f++)
		{
			for(int h = 0; h < (*f)->points.size(); h++)
			{
				double dist = ((*f)->points[h]->pos - corner).norm();

				if(dist < voxel_size * 0.1)
				{
					(*f)->points[h] = cpnts.back();
				}
			}
		}
	}

	// Build voxel graph
	for(int i = 0; i < (int)voxler->voxels.size(); i++)
	{
		StdVector<int> edges = voxler->cornerIndices[i];

		for(int g = 0; g < edges.size(); g++)
		{
			for(int h = g+1; h < edges.size(); h++)
			{
				voxelGraph.AddEdge(edges[g], edges[h], 1);
			}
		}
	}

	// For manual control
	selectedCorner = -1;

	sigmaControl = 1.0;
}

void VoxelDeformer::update()
{
	if(!cpnts.size()) return;

	Vec3d delta = cpnts[selectedCorner]->pos - startPos[selectedCorner];

	StdVector<double> dists = computeNormalizedDistance(selectedCorner);
	
	double sigma = sigmaControl / sqrt(2 * M_PI);

	for(int i = 0; i < cpnts.size(); i++)
	{
		double weight = gaussianFunction(dists[i], 0, sigma);
		cpnts[i]->pos = startPos[i] + (delta * weight);
	}

	foreach(FFD * f, ffd)
	{
		// Compute deformed mesh points
		std::map<int, Point> newPnts = f->applyFixed();

		for(std::map<int, Point>::iterator it = newPnts.begin(); it != newPnts.end(); it++)
		{
			int vid = it->first;
			Surface_mesh::Vertex vert(vid);

			Point pos = it->second;

			mesh->setVertexPos(vert, pos);
		}
	}

	emit(meshDeformed());
}

void VoxelDeformer::push( Vec3d from, Vec3d to )
{
	if(!cpnts.size()) return;

	// Find start

	// Compute weights

	update();
}

void VoxelDeformer::draw()
{
	if(!cpnts.size()) return;

	//voxler->draw();

	// Debug points
	foreach(Point p, debugPoints)	SimpleDraw::IdentifyPoint(p, 1,0,0);
	foreach(Point p, debugPoints2)	SimpleDraw::IdentifyPoint(p, 0,1,0);
	foreach(Point p, debugPoints3)	SimpleDraw::IdentifyPoint(p, 0,0,1);

	// Visualize selected
	if(RANGE(selectedCorner, 0, cpnts.size() - 1))
	{
		StdVector<double> dists = computeNormalizedDistance(selectedCorner);
		double sigma = sigmaControl / sqrt(2 * M_PI);

		uchar rgb[4];
		glPointSize(8);
		glDisable(GL_LIGHTING);
		glBegin(GL_POINTS);
		for(int i = 0; i < cpnts.size(); i++)
		{
			ColorMap::jetColorMap(rgb, gaussianFunction(dists[i], 0, sigma), 0, 1);
			rgb[3] = (1.0 - dists[i]) * 255;
			glColor4ubv(rgb);
			glVertex3dv(cpnts[i]->pos);
		}
		glEnd();
	}

	// Show control points
	StdVector<Point> allPnts;
	foreach(FFD * f, ffd)
		foreach(QControlPoint * cp, f->points)
		allPnts.push_back(cp->pos);

	if(selectedCorner == -1)
		SimpleDraw::IdentifyPoints(allPnts, Color(0.1,0,1,0.5), 4);
	else
		SimpleDraw::IdentifyPoints(allPnts, Color(0.1,0,1,0.1), 2);
}

StdVector<double> VoxelDeformer::computeNormalizedDistance(int startCornerIndex)
{
	StdVector<double> dists;
	double maxDist = DBL_MIN;

	for(int i = 0; i < cpnts.size(); i++)
	{
		std::list<int> path = voxelGraph.DijkstraShortestPath(startCornerIndex, i);
		dists.push_back(path.size());
		maxDist = Max(maxDist, dists.back());
	}

	// Normalize
	for(int i = 0; i < dists.size(); i++)
		dists[i] /= maxDist;

	return dists;
}

void VoxelDeformer::drawNames()
{
	glPointSize(20);

	for(int i = 0; i < cpnts.size(); i++)
	{
		glPushName(i);
		glBegin(GL_POINTS);
		glVertex3dv(cpnts[i]->pos);
		glEnd();
		glPopName();
	}
}

QControlPoint * VoxelDeformer::getQControlPoint( int selected )
{
	if(RANGE(selected, 0, cpnts.size() - 1))
		return cpnts[selected];
	else
		return NULL;
}

void VoxelDeformer::select( int i )
{
	selectedCorner = i;

	startPos.clear();

	foreach(QControlPoint * cp, cpnts)
		startPos.push_back(cp->pos);
}
