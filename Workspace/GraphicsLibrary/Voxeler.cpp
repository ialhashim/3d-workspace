#include "Voxeler.h"
#include "SimpleDraw.h"
#include "Stats.h"

Voxeler::Voxeler( Mesh * src_mesh, double voxel_size )
{
	this->mesh = src_mesh;
	this->voxelSize = voxel_size;

	if(src_mesh == NULL)
		return;

	printf("Computing voxels.."); CreateTimer(timer);

	// For each face in mesh
	for(StdList<Face>::iterator f = mesh->face.begin(); f != mesh->face.end(); f++)
	{
		Face * face = &(*f);

		FaceBounds fb = findFaceBounds( face );

		for(int x = fb.minX; x <= fb.maxX; x++)
		{
			for(int y = fb.minY; y <= fb.maxY; y++)
			{
				for(int z = fb.minZ; z <= fb.maxZ; z++)
				{
					Voxel v(x,y,z);

					if(isVoxelIntersects(v, face) && !kd.has(x,y,z))
					{
						kd.insert3(v.x, v.y, v.z , 1);
						voxels.push_back( v );
					}
				}
			}
		}
	}

	// Compute bounds
	computeBounds();

	// Setup visualization
	setupDraw();

	printf(".voxel count = %d. Done (%d ms).\n", (int)voxels.size(), (int)timer.elapsed());

	// Inner / outer computation
	printf("Computing inside, outside..");
	
	fillInsideOut(innerVoxels, outerVoxels);

	printf("done.");
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

FaceBounds Voxeler::findFaceBounds( Face * f )
{
	FaceBounds fb;

	double minx = 0, miny = 0, minz = 0;
	double maxx = 0, maxy = 0, maxz = 0;

	minx = maxx = f->vec(0).x;
	miny = maxy = f->vec(0).y;
	minz = maxz = f->vec(0).z;

	for(int v = 0; v < 3; v++)
	{
		Vec vec = f->vec(v);

		if (vec.x < minx) minx = vec.x;
		if (vec.x > maxx) maxx = vec.x;

		if (vec.y < miny) miny = vec.y;
		if (vec.y > maxy) maxy = vec.y;

		if (vec.z < minz) minz = vec.z;
		if (vec.z > maxz) maxz = vec.z;
	}

	fb.minX = floor(minx / voxelSize);
	fb.minY = floor(miny / voxelSize);
	fb.minZ = floor(minz / voxelSize);

	fb.maxX = ceil(maxx / voxelSize);
	fb.maxY = ceil(maxy / voxelSize);
	fb.maxZ = ceil(maxz / voxelSize);

	return fb;
}

bool Voxeler::isVoxelIntersects( const Voxel& v, Face * f )
{
	Vec center = Vec(v.x * voxelSize, v.y * voxelSize, v.z * voxelSize);

	double s = voxelSize * 0.5;

	BoundingBox b(center, s,s,s);

	return b.containsTriangle(f->vec(0), f->vec(1), f->vec(2));
}

void Voxeler::draw()
{
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

	Vector<Voxel> temp1, temp2;// = fillOther();

	Vector<double*> insideVoxels = innerVoxels.getAll();
	Vector<double*> outsideVoxels = outerVoxels.getAll();

	foreach(double * pos, insideVoxels)
		temp1.push_back(Voxel(pos[0], pos[1], pos[2]));

	foreach(double * pos, outsideVoxels)
		temp2.push_back(Voxel(pos[0], pos[1], pos[2]));

	for(int i = 0; i < (int) temp1.size(); i++){
		Vec c = temp1[i];
		c *= voxelSize;
		SimpleDraw::DrawSolidBox(c, voxelSize, voxelSize, voxelSize, 0,0,1);
	}

	/*for(int i = 0; i < (int) temp2.size(); i++){
		Vec c = temp2[i];
		c *= voxelSize;
		if(temp2[i].x < 0)
			SimpleDraw::DrawSolidBox(c, voxelSize, voxelSize, voxelSize, 1,0,0);
	}*/

	glEnable(GL_LIGHTING);
}

void Voxeler::setupDraw()
{
	double s = voxelSize * 0.5;
	int n = (int)voxels.size();

	Vector<Vec> c1(n), c2(n), c3(n), c4(n);
	Vector<Vec> bc1(n), bc2(n), bc3(n), bc4(n);

	// Find corners
	for(int i = 0; i < (int)voxels.size(); i++)
	{
		Vec c = voxels[i];	c *= voxelSize;
		c1[i] = Vec(s, s, s) + c; c2[i] = Vec(-s, s, s) + c;
		c3[i] = Vec(-s, -s, s) + c; c4[i] = Vec(s, -s, s) + c;
		bc1[i] = Vec(s, s, -s) + c; bc2[i] = Vec(-s, s, -s) + c;
		bc3[i] = Vec(-s, -s, -s) + c; bc4[i] = Vec(s, -s, -s) + c;
	}

	d1 = glGenLists(1);

	// Faces
	glNewList(d1, GL_COMPILE);
	glBegin(GL_QUADS);
	for(int i = 0; i < (int)voxels.size(); i++)
	{
		// top, left, right, bottom
		gln(0,0,1); glv(c1[i]); glv(c2[i]); glv(c3[i]); glv(c4[i]);
		gln(0,1,0); glv(c1[i]); glv(c2[i]); glv(bc2[i]); glv(bc1[i]);
		gln(0,-1,0); glv(c3[i]); glv(c4[i]); glv(bc4[i]); glv(bc3[i]);
		gln(0,0,-1); glv(bc1[i]); glv(bc2[i]); glv(bc3[i]); glv(bc4[i]);

		// front, back
		gln(1,0,0); glv(c1[i]); glv(c4[i]); glv(bc4[i]); glv(bc1[i]);
		gln(-1,0,0); glv(c2[i]); glv(c3[i]); glv(bc3[i]); glv(bc2[i]);
	}
	glEnd();
	glEndList();

	d2 = glGenLists(1);

	// Lines
	glNewList(d2, GL_COMPILE);
	glBegin(GL_LINES);
	for(int i = 0; i < (int)voxels.size(); i++)
	{
		glv(c1[i]);glv(bc1[i]);glv(c2[i]);glv(bc2[i]);
		glv(c3[i]);glv(bc3[i]);glv(c4[i]);glv(bc4[i]);
		glv(c1[i]);glv(c2[i]);glv(c3[i]);glv(c4[i]);
		glv(c1[i]);glv(c4[i]);glv(c2[i]);glv(c3[i]);
		glv(bc1[i]);glv(bc2[i]);glv(bc3[i]);glv(bc4[i]);
		glv(bc1[i]);glv(bc4[i]);glv(bc2[i]);glv(bc3[i]);
	}
	glEnd();
	glEndList();
}

Vector<Voxel> Voxeler::fillOther()
{
	Vector<Voxel> filled;

	for(int x = minVox.x - 1; x <= maxVox.x + 1; x++)
	{
		for(int y = minVox.y - 1; y <= maxVox.y + 1; y++)
		{
			for(int z = minVox.z - 1; z <= maxVox.z + 1; z++)
			{
				if(!kd.has(x,y,z))
				{
					filled.push_back(Voxel(x,y,z));
				}
			}
		}
	}

	return filled;
}

void Voxeler::fillInsideOut(KDTree & inside, KDTree & outside)
{
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
	Stack<Voxel> stack;

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
