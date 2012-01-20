#include "Voxeler.h"
#include "SimpleDraw.h"
#include "Stats.h"

Voxeler::Voxeler( QSurfaceMesh * src_mesh, double voxel_size, bool verbose /*= false*/ )
{
	this->mesh = src_mesh;
	this->voxelSize = voxel_size;
	this->isVerbose = verbose;
	this->isReadyDraw = false;

	if(mesh == NULL)
		return;

	if(isVerbose) printf("Computing voxels..");

	mesh->assignFaceArray();

	// For each face in mesh
	foreach(Surface_mesh::Face f, mesh->face_array)
	{
		FaceBounds fb = findFaceBounds( f );

		for(int x = fb.minX; x <= fb.maxX; x++)
		{
			for(int y = fb.minY; y <= fb.maxY; y++)
			{
				for(int z = fb.minZ; z <= fb.maxZ; z++)
				{
					Voxel v(x,y,z);

					if(isVoxelIntersects(v, f) && !kd.has(x,y,z)){
						kd.insert3(v.x, v.y, v.z, voxels.size());
						voxels.push_back( v );
					}
				}
			}
		}
	}
	
	if(isVerbose) printf(".voxel count = %d.\n", (int)voxels.size());

	// Inner / outer computation
	//fillInsideOut(innerVoxels, outerVoxels);
	
	if(isVerbose) printf("done.");
}

void Voxeler::update()
{
	// Compute bounds
	computeBounds();

	// Setup visualization
	setupDraw();
}

void Voxeler::computeBounds()
{
	minVox = Voxel(INT_MAX, INT_MAX, INT_MAX);
	maxVox = Voxel(INT_MIN, INT_MIN, INT_MIN);

	for(int i = 0; i < (int)voxels.size(); i++)
	{
		Voxel v = voxels[i];

		minVox.toMin(v);
		maxVox.toMax(v);
	}
}

FaceBounds Voxeler::findFaceBounds( QSurfaceMesh::Face f )
{
	FaceBounds fb;

	double minx = 0, miny = 0, minz = 0;
	double maxx = 0, maxy = 0, maxz = 0;

	std::vector<Vec3d> f_vec =  mesh->facePoints(f);

	minx = maxx = f_vec[0].x();
	miny = maxy = f_vec[0].y();
	minz = maxz = f_vec[0].z();

	for(int v = 0; v < 3; v++)
	{
		Vec3d vec = f_vec[v];

		if (vec.x() < minx) minx = vec.x();
		if (vec.x() > maxx) maxx = vec.x();

		if (vec.y() < miny) miny = vec.y();
		if (vec.y() > maxy) maxy = vec.y();

		if (vec.z() < minz) minz = vec.z();
		if (vec.z() > maxz) maxz = vec.z();
	}

	fb.minX = floor(minx / voxelSize);
	fb.minY = floor(miny / voxelSize);
	fb.minZ = floor(minz / voxelSize);

	fb.maxX = ceil(maxx / voxelSize);
	fb.maxY = ceil(maxy / voxelSize);
	fb.maxZ = ceil(maxz / voxelSize);

	return fb;
}

bool Voxeler::isVoxelIntersects( const Voxel& v, QSurfaceMesh::Face f )
{
	Vec3d center = Vec3d(v.x * voxelSize, v.y * voxelSize, v.z * voxelSize);

	double s = voxelSize * 0.5;

	BoundingBox b(center, s,s,s);

	std::vector<Vec3d> f_vec =  mesh->facePoints(f);

	return b.containsTriangle(f_vec[0], f_vec[1], f_vec[2]);
}

void Voxeler::draw()
{
	if(!isReadyDraw)
		update();

	glEnable(GL_LIGHTING);
	glShadeModel(GL_FLAT);

	glColor3f(0, 0.8f, 0);
	//glCallList(d1);

	glDisable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);

	glColor3f(0, 1, 0);
	glLineWidth(1.0);
	glCallList(d2);

	// DEBUG == DELETE ME:

	std::vector<Voxel> temp1, temp2;// = fillOther();

	std::vector<double*> insideVoxels = innerVoxels.getAll();
	std::vector<double*> outsideVoxels = outerVoxels.getAll();

	foreach(double * pos, insideVoxels)
		temp1.push_back(Voxel(pos[0], pos[1], pos[2]));

	foreach(double * pos, outsideVoxels)
		temp2.push_back(Voxel(pos[0], pos[1], pos[2]));

	for(int i = 0; i < (int) temp1.size(); i++){
		Vec3d c = temp1[i];
		c *= voxelSize;
		SimpleDraw::DrawSolidBox(c, voxelSize, voxelSize, voxelSize, 0,0,1);
	}

	/*for(int i = 0; i < (int) temp2.size(); i++){
		Vec3d c = temp2[i];
		c *= voxelSize;
		if(temp2[i].x() < 0)
			SimpleDraw::DrawSolidBox(c, voxelSize, voxelSize, voxelSize, 1,0,0);
	}*/

	glEnable(GL_LIGHTING);
}

void Voxeler::setupDraw()
{
	double s = voxelSize * 0.5;
	int n = (int)voxels.size();
	
	std::vector< std::vector<Vec3d> > corner(n, std::vector<Vec3d>(8));

	// Find corners
	for(int i = 0; i < n; i++)
	{
		Vec3d c = voxels[i];	c *= voxelSize;
		corner[i][0] = Vec3d(s, s, s) + c;		corner[i][1] = Vec3d(-s, s, s) + c;
		corner[i][2] = Vec3d(-s, -s, s) + c;	corner[i][3] = Vec3d(s, -s, s) + c;
		corner[i][4] = Vec3d(s, s, -s) + c;		corner[i][5] = Vec3d(-s, s, -s) + c;
		corner[i][6] = Vec3d(-s, -s, -s) + c;	corner[i][7] = Vec3d(s, -s, -s) + c;
	}

	// Save corners
	corners.clear();
	cornerIndices.clear();
	cornerCorrespond.resize(n*8);

	for(int i = 0; i < n; i++)
	{
		std::vector<int> p(8, 0);

		for(int j = 0; j < 8; j++)
		{
			if(!corner_kd.has(corner[i][j])){
				p[j] = corners.size();
				corners.push_back(corner[i][j]);
				corner_kd.insert3(corner[i][j][0], corner[i][j][1], corner[i][j][2], p[j]);
			}else p[j] = corner_kd.getData(&corner[i][j][0]);

			cornerCorrespond[p[j]] = i; // will belong to last one
		}

		cornerIndices.push_back(p);
	}

	d1 = glGenLists(1);

	// Faces
	glNewList(d1, GL_COMPILE);
	glBegin(GL_QUADS);
	for(int i = 0; i < (int)voxels.size(); i++)
	{
		// top, left, right, bottom
		gln(0,0,1); glv(corner[i][0]); glv(corner[i][1]); glv(corner[i][2]); glv(corner[i][3]);
		gln(0,1,0); glv(corner[i][0]); glv(corner[i][1]); glv(corner[i][5]); glv(corner[i][4]);
		gln(0,-1,0); glv(corner[i][2]); glv(corner[i][3]); glv(corner[i][7]); glv(corner[i][6]);
		gln(0,0,-1); glv(corner[i][4]); glv(corner[i][5]); glv(corner[i][6]); glv(corner[i][7]);

		// front, back
		gln(1,0,0); glv(corner[i][0]); glv(corner[i][3]); glv(corner[i][7]); glv(corner[i][4]);
		gln(-1,0,0); glv(corner[i][1]); glv(corner[i][2]); glv(corner[i][6]); glv(corner[i][5]);
	}
	glEnd();
	glEndList();

	d2 = glGenLists(1);

	// Lines
	glNewList(d2, GL_COMPILE);
	glBegin(GL_LINES);
	for(int i = 0; i < (int)voxels.size(); i++)
	{
		glv(corner[i][0]);glv(corner[i][4]);glv(corner[i][1]);glv(corner[i][5]);
		glv(corner[i][2]);glv(corner[i][6]);glv(corner[i][3]);glv(corner[i][7]);
		glv(corner[i][0]);glv(corner[i][1]);glv(corner[i][2]);glv(corner[i][3]);
		glv(corner[i][0]);glv(corner[i][3]);glv(corner[i][1]);glv(corner[i][2]);
		glv(corner[i][4]);glv(corner[i][5]);glv(corner[i][6]);glv(corner[i][7]);
		glv(corner[i][4]);glv(corner[i][7]);glv(corner[i][5]);glv(corner[i][6]);
	}
	glEnd();
	glEndList();

	isReadyDraw = true;
}

void Voxeler::drawVoxels( const std::vector< Voxel > & voxels, double voxel_size )
{
	double s = voxel_size * 0.5;
	int n = (int)voxels.size();

	std::vector<Vec3d> c1(n), c2(n), c3(n), c4(n);
	std::vector<Vec3d> bc1(n), bc2(n), bc3(n), bc4(n);

	// Find corners
	for(int i = 0; i < (int)voxels.size(); i++){
		Vec3d c = voxels[i];	c *= voxel_size;
		c1[i] = Vec3d(s, s, s) + c; c2[i] = Vec3d(-s, s, s) + c;
		c3[i] = Vec3d(-s, -s, s) + c; c4[i] = Vec3d(s, -s, s) + c;
		bc1[i] = Vec3d(s, s, -s) + c; bc2[i] = Vec3d(-s, s, -s) + c;
		bc3[i] = Vec3d(-s, -s, -s) + c; bc4[i] = Vec3d(s, -s, -s) + c;
	}

	glColor3d(1,0,0);
	glLineWidth(3.0f);
	glDisable(GL_LIGHTING);

	// Lines
	glBegin(GL_LINES);
	for(int i = 0; i < (int)voxels.size(); i++){
		glv(c1[i]);glv(bc1[i]);glv(c2[i]);glv(bc2[i]);
		glv(c3[i]);glv(bc3[i]);glv(c4[i]);glv(bc4[i]);
		glv(c1[i]);glv(c2[i]);glv(c3[i]);glv(c4[i]);
		glv(c1[i]);glv(c4[i]);glv(c2[i]);glv(c3[i]);
		glv(bc1[i]);glv(bc2[i]);glv(bc3[i]);glv(bc4[i]);
		glv(bc1[i]);glv(bc4[i]);glv(bc2[i]);glv(bc3[i]);
	}
	glEnd();

	glEnable(GL_LIGHTING);
}

std::vector<Voxel> Voxeler::fillOther()
{
	std::vector<Voxel> filled;

	for(int x = minVox.x - 1; x <= maxVox.x + 1; x++){
		for(int y = minVox.y - 1; y <= maxVox.y + 1; y++){
			for(int z = minVox.z - 1; z <= maxVox.z + 1; z++){
				if(!kd.has(x,y,z))
					filled.push_back(Voxel(x,y,z));
			}
		}
	}

	return filled;
}

void Voxeler::fillInsideOut(KDTree & inside, KDTree & outside)
{
	printf("Computing inside, outside..");

	fillOuter(outside);

	// Compute inner as complement of outside
	for(int x = minVox.x - 1; x <= maxVox.x + 1; x++){
		for(int y = minVox.y - 1; y <= maxVox.y + 1; y++){
			for(int z = minVox.z - 1; z <= maxVox.z + 1; z++){
				if(!kd.has(x,y,z) && !outside.has(x,y,z)){
					inside.insert3(x,y,z,1);
				}
			}
		}
	}
}

void Voxeler::fillOuter(KDTree & outside)
{
	std::stack<Voxel> stack;

	stack.push(maxVox + Voxel(1,1,1));

	while(!stack.empty())
	{
		// Get next square
		Voxel c = stack.top(); // get current voxel
		stack.pop();

		// Base case:
		if( !kd.has(c.x, c.y, c.z) && !outside.has(c.x, c.y, c.z) )
		{
			// Otherwise, add it to set of outside voxels
			outside.insert3(c.x, c.y, c.z, 1);

			// Visit neighbors
			if(c.x < maxVox.x + 1) stack.push( c + Voxel( 1, 0, 0) );
			if(c.y < maxVox.y + 1) stack.push( c + Voxel( 0, 1, 0) );
			if(c.z < maxVox.z + 1) stack.push( c + Voxel( 0, 0, 1) );

			if(c.x > minVox.x - 1) stack.push( c + Voxel(-1, 0, 0) );
			if(c.y > minVox.y - 1) stack.push( c + Voxel( 0,-1, 0) );
			if(c.z > minVox.z - 1) stack.push( c + Voxel( 0, 0,-1) );
		}
	}
}

std::vector<Voxel> Voxeler::Intersects(Voxeler * other)
{
	std::vector<Voxel> intersection;

	Voxeler *minVoxeler = this, *maxVoxeler = other;

	// Swap with minimum
	if(other->voxels.size() < this->voxels.size()){
		minVoxeler = other;
		maxVoxeler = this;
	}

	for(int i = 0; i < minVoxeler->voxels.size(); i++)
	{
		Voxel v = minVoxeler->voxels[i];

		if(maxVoxeler->kd.has(v.x, v.y, v.z))
			intersection.push_back(v);
	}

	return intersection;
}

std::map<int, Voxel> Voxeler::around(Point p)
{
	std::map<int, Voxel> result;

	int x = p.x() / voxelSize;
	int y = p.y() / voxelSize;
	int z = p.z() / voxelSize;

	for(int i = -1; i <= 1; i += 1){
		for(int j = -1; j <= 1; j += 1){
			for(int k = -1; k <= 1; k += 1){
				Voxel v(x + i, y + j, z + k);

				Vec3d vpos(v.x, v.y, v.z);

				if(kd.has(v.x, v.y, v.z)){
					result[kd.getData(&vpos[0])] = v;
				}
			}
		}
	}

	return result;
}

void Voxeler::grow()
{
	int N = (int)voxels.size();

	for(int i = 0; i < N; i++)
	{
		Voxel curVoxel = voxels[i];

		for(int x = -1; x <= 1; x += 1){
			for(int y = -1; y <= 1; y++){
				for(int z = -1; z <= 1; z++){
					Voxel v(curVoxel.x + x, curVoxel.y + y, curVoxel.z + z);

					if(!kd.has(v.x,v.y,v.z))
					{
						kd.insert3(v.x, v.y, v.z, voxels.size());
						voxels.push_back( v );
					}
				}
			}
		}
	}

	printf("\nVoxler grown from (%d) to (%d).\n", N, (int)voxels.size());

	isReadyDraw = false;
}

std::vector< Point > Voxeler::getCorners( int vid )
{
	std::vector< Point > result;

	for(int i = 0; i < 8; i++)
		result.push_back(corners[cornerIndices[vid][i]]);

	return result;
}

int Voxeler::getClosestVoxel( Vec3d point )
{
	return cornerCorrespond[ corner_kd.getData(&point[0]) ];
}

int Voxeler::getEnclosingVoxel( Vec3d point )
{
	int N = (int)voxels.size();

	double s = voxelSize * 0.5;

	double minDist = DBL_MAX;
	int closestVoxel = -1;

	for(int i = 0; i < N; i++)
	{
		Voxel curVoxel = voxels[i];

		Point voxelCenter(curVoxel.x * s, curVoxel.y * s, curVoxel.z * s);

		double curDist = (voxelCenter - point).norm();

		if(curDist < minDist){
			closestVoxel = i;
			minDist = curDist;
		}
	}

	return closestVoxel;
}

int Voxeler::getVoxelIndex( Voxel v )
{
	for(int i = 0; i < voxels.size(); i++)
	{
		Voxel w = voxels[i];

		if(w.x == v.x && w.y == v.y && w.z == v.z)
			return i;
	}

	return -1;
}
