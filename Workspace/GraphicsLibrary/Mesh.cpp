#include "Mesh.h"

#include "SimpleDraw.h"		// Debug helper

#include <omp.h>			// OpenMP

#include "Graph.h"

Mesh::Mesh(int expectedNumVerts)
{
	
	// Defaults
	isVisible = true;
	isDrawSmooth = true;

	isReady = false;
	isDrawWireframe = false;
	isDrawVertices = false;
	isTransparent = false;
	isFlatShade = false;
	isDrawAsPoints = false;
	isShowVertexNormals = false;
	isShowFaceNormals = false;

	this->vbo = NULL;
	this->octree = NULL;

	this->radius = 0.0;
	this->normalize_scale = 1.0;

	this->vertex.reserve(expectedNumVerts);
	this->vertexInfo.reserve(expectedNumVerts);
	this->vNormal.reserve(expectedNumVerts);

	//printf("Empty mesh created.\n\n");

	// debug items
	selectedFace = -1;
	selectedVertex = -1;
}

Mesh::Mesh(const Mesh& fromMesh)
{
	this->octree = NULL;
	this->center = fromMesh.center;

	this->id = fromMesh.id;
	this->isReady = fromMesh.isReady;
	this->maxBound = fromMesh.maxBound;
	this->minBound = fromMesh.minBound;
	this->radius = fromMesh.radius;
	this->normalize_scale = fromMesh.normalize_scale;

	this->face = fromMesh.face;
	this->fNormal = fromMesh.fNormal;

	this->vertex = fromMesh.vertex;
	this->vertexInfo = Vector<VertexDetail>(fromMesh.vertexInfo.size());
	this->vNormal = fromMesh.vNormal;
	this->vColor = fromMesh.vColor;

	int v1, v2, v3;

	for(StdList<Face>::iterator f = this->face.begin(); f != this->face.end(); f++)
	{
		v1 = f->vIndex[0];
		v2 = f->vIndex[1];
		v3 = f->vIndex[2];

		f->v[0] = &this->vertex[v1];
		f->v[1] = &this->vertex[v2];
		f->v[2] = &this->vertex[v3];

		this->vertexInfo[v1].insertFace(&(*f)); 
		this->vertexInfo[v1].index = v1;

		this->vertexInfo[v2].insertFace(&(*f));
		this->vertexInfo[v2].index = v2;

		this->vertexInfo[v3].insertFace(&(*f));
		this->vertexInfo[v3].index = v3;

		this->faceIndexMap[f->index] = &(*f);
	}

	if(fromMesh.vbo == NULL)
		this->vbo = NULL;
	else
		this->vbo = new VBO(&this->vertex, &this->vNormal, &this->vColor, &this->face);

	//this->tempUmbrellas = fromMesh.tempUmbrellas;
	this->tempUmbrellas.clear();

	this->isReady = fromMesh.isReady;
	this->isVisible = fromMesh.isVisible;
	this->isDrawSmooth = fromMesh.isDrawSmooth;
	this->isDrawWireframe = fromMesh.isDrawWireframe;
	this->isDrawVertices = fromMesh.isDrawVertices;
	this->isTransparent = fromMesh.isTransparent;
	this->isFlatShade = fromMesh.isFlatShade;
	this->isDrawAsPoints = fromMesh.isDrawAsPoints;
	this->isShowVertexNormals = fromMesh.isShowVertexNormals;
	this->isShowFaceNormals = fromMesh.isShowFaceNormals;

	// debug items
	selectedFace = -1;
	selectedVertex = -1;
}

Mesh::~Mesh()
{
	this->vertex.clear();
	this->vertexInfo.clear();
	this->face.clear();

	this->vNormal.clear();
	this->vColor.clear();
	this->fNormal.clear();

	if(octree)
		delete octree;
	if(vbo)
		delete vbo;
}

Mesh& Mesh::operator= (const Mesh& fromMesh)
{
	if (this != &fromMesh) {
		this->octree = fromMesh.octree;

		this->center = fromMesh.center;

		this->id = fromMesh.id;
		this->isReady = fromMesh.isReady;
		this->maxBound = fromMesh.maxBound;
		this->minBound = fromMesh.minBound;
		this->radius = fromMesh.radius;
		this->normalize_scale = fromMesh.normalize_scale;

		this->face = fromMesh.face;
		this->fNormal = fromMesh.fNormal;
		this->faceIndexMap = fromMesh.faceIndexMap;

		this->vertex = fromMesh.vertex;
		this->vertexInfo = Vector<VertexDetail>(fromMesh.vertexInfo.size());
		this->vNormal = fromMesh.vNormal;
		this->vColor = fromMesh.vColor;

		int v1, v2, v3;

		for(StdList<Face>::iterator f = this->face.begin(); f != this->face.end(); f++)
		{
			v1 = f->vIndex[0];
			v2 = f->vIndex[1];
			v3 = f->vIndex[2];

			f->v[0] = &this->vertex[v1];
			f->v[1] = &this->vertex[v2];
			f->v[2] = &this->vertex[v3];

			this->vertexInfo[v1].insertFace(&(*f)); 
			this->vertexInfo[v1].index = v1;

			this->vertexInfo[v2].insertFace(&(*f));
			this->vertexInfo[v2].index = v2;

			this->vertexInfo[v3].insertFace(&(*f));
			this->vertexInfo[v3].index = v3;

			this->faceIndexMap[f->index] = &(*f);
		}

		if(fromMesh.vbo == NULL)
			this->vbo = NULL;
		else
			this->vbo = new VBO(&this->vertex, &this->vNormal, &this->vColor, &this->face);

		//this->tempUmbrellas = fromMesh.tempUmbrellas;
		this->tempUmbrellas.clear();

		this->isReady = fromMesh.isReady;
		this->isVisible = fromMesh.isVisible;
		this->isDrawSmooth = fromMesh.isDrawSmooth;
		this->isDrawWireframe = fromMesh.isDrawWireframe;
		this->isDrawVertices = fromMesh.isDrawVertices;
		this->isTransparent = fromMesh.isTransparent;
		this->isFlatShade = fromMesh.isFlatShade;
		this->isDrawAsPoints = fromMesh.isDrawAsPoints;
		this->isShowVertexNormals = fromMesh.isShowVertexNormals;
		this->isShowFaceNormals = fromMesh.isShowFaceNormals;

		// debug items
		selectedFace = -1;
		selectedVertex = -1;
	}

	return *this;
}

void Mesh::mergeWith(const Mesh& other)
{
	this->isReady = false;

	int offset = this->numberOfVertices();

	for(Vector<Point3D>::const_iterator it = other.vertex.begin(); it != other.vertex.end(); it++)
	{
		addVertex(*it, this->vertex.size());

		vColor.push_back(Color4());
	}

	for(StdList<Face>::const_iterator it = other.face.begin(); it != other.face.end(); it++)
	{
		addFace(it->vIndex[0] + offset, it->vIndex[1] + offset, it->vIndex[2] + offset, this->face.size());
	}
}

void Mesh::addVertex(double x, double y, double z, int index)
{
	vertex.push_back(Point3D(x, y, z));
	vertexInfo.push_back(VertexDetail(index));
}

void Mesh::addVertex(Vec v, int index)
{
	addVertex(v.x, v.y, v.z, index);
}

void Mesh::addFace(int v0, int v1, int v2, int index, bool forceOrientation)
{
	if(forceOrientation)
	{
		// if vertex normals are known, we can find correct orientation
		Vec newNormal = (vertex[v1] - vertex[v0]) ^ (vertex[v2] - vertex[v0]);

		// If they are not pointing at the same direction, swap first and last
		if(newNormal * computeVNormalAt(v0) < 0)
		{
			int temp = v1;
			v1 = v2;
			v2 = temp;
		}
	}

	face.push_back(Face(v0, v1, v2, &this->vertex[v0], &this->vertex[v1], &this->vertex[v2], index));

	if(v0 == v1) printf("WARNING: degenerate face (v1,v2) \n");
	if(v0 == v2) printf("WARNING: degenerate face (v1,v3) \n");
	if(v1 == v2) printf("WARNING: degenerate face (v2,v3) \n");

	Face * f = &face.back();

	faceIndexMap[index] = f;

	vertexInfo[v0].insertFace(f);
	vertexInfo[v1].insertFace(f);
	vertexInfo[v2].insertFace(f);
}

void Mesh::computeBounds()
{
	double minx, miny, minz, maxx, maxy, maxz;

	minx = maxx = vertex[0].x;
	miny = maxy = vertex[0].y;
	minz = maxz = vertex[0].z;

	int N = this->numberOfVertices();

	for(int i=0; i < N; i++)
	{
		if (vertex[i].x < minx) minx = vertex[i].x;
		if (vertex[i].x > maxx) maxx = vertex[i].x;
		if (vertex[i].y < miny) miny = vertex[i].y;
		if (vertex[i].y > maxy) maxy = vertex[i].y;
		if (vertex[i].z < minz) minz = vertex[i].z;
		if (vertex[i].z > maxz) maxz = vertex[i].z;
	}

	minBound = Vec(minx, miny, minz);
	maxBound = Vec(maxx, maxy, maxz);

	center = (minBound + maxBound) / 2.0;
	radius = (maxBound - center).norm();
}

void Mesh::moveToCenter()
{
	computeBounds();

	int N = this->numberOfVertices();

	for(int i=0; i < N; i++)
		vertex[i] -= center;

	computeBounds();
}

void Mesh::normalizeScale()
{
	double sizeX = maxBound[0] - minBound[0];
	double sizeY = maxBound[1] - minBound[1];
	double sizeZ = maxBound[2] - minBound[2];

	maxScale = Max(sizeX, Max(sizeY, sizeZ));

	double scale = normalize_scale / maxScale;

	int N = this->numberOfVertices();

#pragma omp parallel for
	for(int i=0; i < N; i++)
	{
		vertex[i].x = (vertex[i].x - center.x) * scale;
		vertex[i].y = (vertex[i].y - center.y) * scale;
		vertex[i].z = (vertex[i].z - center.z) * scale;
	}

	computeBounds();

	printf("Mesh Radius = %f (scaled:%f)\n", radius, scale);

	setDirtyVBO(true);
}

void Mesh::computeNormals()
{
	// Compute face normals
	fNormal = Vector<Normal>(face.size());

	for(StdList<Face>::iterator f = this->face.begin(); f != this->face.end(); f++)
		fNormal[f->index].set(f->normal());

	int N = vertex.size();

	// Compute vertex normals
	vNormal = Vector<Normal>(N);

	for( int i = 0; i < N; ++i )
	{
		Vector<Face *> currFace = vertexInfo[i].ifaces;
		int numFaces = currFace.size();

		if(numFaces > 0)
		{
			Normal v_normal;

			for (int j = 0; j < numFaces; j++)
				v_normal += fNormal[currFace[j]->index];

			vNormal[i] = v_normal.unit();
		}

		if(i % (N/5) == 1)	printf(".");
	}
}

Vec Mesh::computeVNormalAt(int vi)
{
	Normal n;

	foreach (Face * f, vertexInfo[vi].ifaces)
		n += f->normal();

	return n.unit();
}

double Mesh::computeVolume()
{
	double vol = 0;

	Point3D p1, p2, p3;
	Point3D g, n;

	for(StdList<Face>::iterator f = this->face.begin(); f != this->face.end(); f++)
		vol += f->volume();

	return vol;
}

double Mesh::computeRadius()
{
	computeBounds();
	return this->radius;
}

Vec Mesh::computeCenter()
{
	Vec centers;

	for(StdList<Face>::iterator f = this->face.begin(); f != this->face.end(); f++)
		centers += f->center();

	return centers / face.size();
}

void Mesh::loadFromFile( const char* fileName )
{
	StdString fileNameString = fileName;
	int n = fileNameString.length();

	StdString fileExt = fileNameString.substr(n - 3, 3);

	if(fileExt == "obj")
		loadFromFileOBJ(fileName);
	else if (fileExt == "off")
		loadFromFileOFF(fileName);
}

void Mesh::loadFromFileOBJ( const char* fileName )
{

	isReady = false;
	opendFileName = fileName;

	CreateTimer(allStartTime);

	printf("Loading...(%s)\n", fileName);

	std::string inputLine;
	FileStream file (fileName);

	float x,y,z;
	int v1,v2,v3;

	int vCount = 0;
	int fCount = 0;

	// Does this help? to reserve 100K vertices
	vertex.reserve(100000);
	vertexInfo.reserve(100000);

	if (file.is_open())
	{
		while (!file.eof())
		{
			GetLine (file, inputLine);

			switch(inputLine[0])
			{
			case 'v':
				if(sscanf(inputLine.c_str(), "v %f %f %f", &x,&y,&z) == 3)
				{
					this->addVertex (x, y, z, vCount++);

					if(vCount % 10000 == 0)	printf(".");
				}
				break;

			case 'f':
				if(sscanf(inputLine.c_str(), "f %d %d %d", &v1,&v2,&v3) == 3 ||
					sscanf(inputLine.c_str(), "f %d//%d %d//%d %d//%d", &v1,&v1 , &v2,&v2 , &v3,&v3) == 6)
				{
					this->addFace(v1 - 1, v2 - 1, v3 - 1, fCount);

					if(fCount % 10000 == 0)	printf(".");
					fCount++;
				}

				break;
			}
		}

		// default color
		vColor = Vector<Color4>(vertex.size(), Color4());

		// we are done, close the file
		file.close();
	}

	// Center mesh into world
	moveToCenter();

	// Normalize mesh
	//normalizeScale();

	// Normals
	printf("\nComputing Normals.."); CreateTimer(normalTimer);
	computeNormals();
	printf("Done (%d ms).\n", (int)normalTimer.elapsed());

	// VBO representation
	printf("Creating VBO object.."); CreateTimer(vboTimer);
	createVBO();
	printf("Done (%d ms).\n", (int)vboTimer.elapsed());

	// Umbrella structure
	printf("Collecting Umbrellas.."); CreateTimer(umbrellasTimer);
	getUmbrellas();
	printf("Done (%d ms).\n", (int)umbrellasTimer.elapsed());

	printf("\n\t V = \t%d\tF = \t%d\n", (int)vertex.size(), (int)face.size());
	printf("\nMesh file loaded successfully. (%d ms)\n", (int)allStartTime.elapsed());

	isReady = true;
}

void Mesh::loadFromFileOFF( const char* fileName )
{
	isReady = false;
	opendFileName = fileName;

	CreateTimer(allStartTime);

	printf("Loading...(%s)\n", fileName);

	std::string inputLine;
	FileStream file (fileName);

	float x,y,z;
	int v1,v2,v3,nv;

	int vCount, fCount, eCount, total_vCount, total_fCount, total_eCount;
	vCount = fCount = eCount = total_vCount = total_fCount = total_eCount = 0;

	bool readingVertices, readingFaces;
	readingVertices = readingFaces = false;

	if (file.is_open())
	{
		while (!file.eof())
		{
			GetLine (file, inputLine);

			// Skip comments, empty lines
			if(inputLine[0] == '#')	
				continue;

			if(readingFaces)
			{
				if(sscanf(inputLine.c_str(), "%d %d %d %d",&nv,&v1,&v2,&v3) == 4)
				{
					this->addFace(v1, v2, v3, fCount++);

					if(fCount % 10000 == 0)	printf(".");

					if(fCount == total_fCount)
						break;
				}
			}

			if(readingVertices && vCount < total_vCount && inputLine.length() > 4)
			{
				if(sscanf(inputLine.c_str(), "%f %f %f",&x,&y,&z) == 3)
				{
					this->addVertex (x, y, z, vCount++);
					if(vCount % 10000 == 0)	printf(".");

					if(vCount == total_vCount)
					{
						readingVertices = false;
						readingFaces = true;
						printf(".Reading faces ");
					}
					continue;
				}
			}

			// Try to find first line with three integers
			if(!readingVertices && !readingFaces && 
				sscanf(inputLine.c_str(), "%d %d %d",&total_vCount,&total_fCount,&total_eCount) == 3)
			{
				readingVertices = true;
				readingFaces = false;
				printf(".Reading vertices ");
			}
		}

		// default color
		vColor = Vector<Color4>(vertex.size(), Color4());

		// we are done, close the file
		file.close();
	}

	// Center mesh into world
	moveToCenter();

	// Normalize mesh
	//normalizeScale();

	// Normals
	printf("\nComputing Normals.."); CreateTimer(normalTimer);
	computeNormals();
	printf("Done (%d ms).\n", (int)normalTimer.elapsed());

	// VBO representation
	printf("Creating VBO object.."); CreateTimer(vboTimer);
	createVBO();
	printf("Done (%d ms).\n", (int)vboTimer.elapsed());

	// Umbrella structure
	printf("Collecting Umbrellas.."); CreateTimer(umbrellasTimer);
	getUmbrellas();
	printf("Done (%d ms).\n", (int)umbrellasTimer.elapsed());

	printf("\n\t V = \t%d\tF = \t%d\n", (int)vertex.size(), (int)face.size());
	printf("\nMesh file loaded successfully. (%d ms)\n", (int)allStartTime.elapsed());

	isReady = true;
}

void Mesh::saveToFile(const char* fileName)
{
	StdString saveFileName;

	if(fileName == NULL)
		saveFileName = id + ".obj";
	else
		saveFileName = fileName;

	FILE *fp = fopen(saveFileName.c_str(),"w");

	if(fp == NULL)	return;

	int numVerts = numberOfVertices();
	int numFaces = numberOfFaces();

	fprintf(fp, "# Vertices %d, Faces %d\n", numVerts, numFaces);

	// Write vertices
	for(int i = 0; i < numVerts; i++){
		Point3D * v = this->v(i);
		fprintf(fp, "v %f %f %f\n", v->x, v->y, v->z);
	}

	// Write triangles
	for(int i = 0; i < numFaces; i++){
		Face * f = this->f(i);

		if(f->flag != FF_INVALID_VINDEX && f->VIndex(0) != -1)
		{
			fprintf(fp, "f %d %d %d\n", f->vIndex[0] + 1, f->vIndex[1] + 1, f->vIndex[2] + 1);
		}
	}

	fclose(fp);
}

void Mesh::createVBO()
{
	this->vbo = new VBO(&vertex, &vNormal, &vColor, &face);
}

void Mesh::updateVBO()
{
	this->vbo->update();
}

int Mesh::numberOfVertices()
{
	return vertex.size();
}

int Mesh::numberOfFaces()
{
	return face.size();
}

bool Mesh::withID(StdString ID)
{
	return this->id == ID;
}

void Mesh::setColor( const Color4 & color )
{
	if(color.a() < 1.0)
		isTransparent = true;
	else
		isTransparent = false;

	vColor = Vector<Color4>(vertex.size(), color);

	if(vbo)	vbo->setDirty(true);
}

void Mesh::setColor(int r, int g, int b, int a)
{
	if(a < 255)
		isTransparent = true;
	else
		isTransparent = false;

	vColor = Vector<Color4>(vertex.size(), Color4(r,g,b,a));

	if(vbo)	vbo->setDirty(true);
}

void Mesh::setColor(Vector<int> & vertices, int r, int g, int b, int a)
{
	if(a < 255)
		isTransparent = true;
	else
		isTransparent = false;

	for(int i = 0; i < (int)vertices.size(); i++)
		vColor[vertices[i]].set(r, g, b, a);

	if(vbo)	vbo->setDirty(true);
}

void Mesh::setColorFace(int faceIndex, int r, int g, int b, int a)
{
	if(a < 255)
		isTransparent = true;
	else
		isTransparent = false;

	Face * f = faceIndexMap[faceIndex];

	vColor[f->VIndex(0)].set(r, g, b, a);
	vColor[f->VIndex(1)].set(r, g, b, a);
	vColor[f->VIndex(2)].set(r, g, b, a);

	if(vbo)	vbo->setDirty(true);
}

void Mesh::setColorFaces(Vector<int> & faces, int r, int g, int b, int a)
{
	if(a < 255)
		isTransparent = true;
	else
		isTransparent = false;

	for(int i = 0; i < (int)faces.size(); i++)
	{
		Face * f = this->f(faces[i]);

		vColor[f->VIndex(0)].set(r, g, b, a);
		vColor[f->VIndex(1)].set(r, g, b, a);
		vColor[f->VIndex(2)].set(r, g, b, a);
	}

	if(vbo) vbo->setDirty(true);
}

void Mesh::colorSplit(Vec & pos, Vec & dir)
{
	Vector<Face *> partA;

	for(StdList<Face>::iterator f = face.begin(); f != face.end(); f++)
	{
		if(f->P(0)->above(pos, dir) > 0)
		{
			vColor[f->VIndex(0)].set(255, 255, 0);
			vColor[f->VIndex(1)].set(255, 255, 0);
			vColor[f->VIndex(2)].set(255, 255, 0);
		}
	}

	if(vbo)	vbo->setDirty(true);
}

void Mesh::clearColors()
{
	vColor = Vector<Color4>(vertex.size(), Color4());
}

void Mesh::colorMarkBorders(int r, int g, int b, int a)
{
	StdMap<int,int> facesMap;
	Vector<int> facesList;

	flagBorderVertices();

	for(Vector<VertexDetail>::iterator vd = vertexInfo.begin(); vd != vertexInfo.end(); vd++)
	{
		if(vd->isBorderFlag())
		{
			for(Vector<Face *>::iterator it = vd->ifaces.begin(); it != vd->ifaces.end(); it++)
				facesMap[(*it)->index] = 1;
		}
	}

	for(StdMap<int,int>::iterator it = facesMap.begin(); it != facesMap.end(); it++)
		facesList.push_back(it->first);

	this->setColorFaces(facesList, r, g, b, a);
}

void Mesh::flagBorderVertices()
{
	clearAllVertexFlag();

	int numBorder = 0;
	int N = vertexInfo.size();

	for(int vi = 0; vi < N; vi++)
	{
		if(vertexInfo[vi].checkIsBorder())
		{
			vertexInfo[vi].flag = VF_BORDER;
			numBorder++;
		}
	}
	//printf(" # at border (%d) ", numBorder);
}

Vector<int> Mesh::getBorderVertices()
{
	clearAllVertexFlag();

	Vector<int> borderVerts;
	borderVerts.reserve(vertex.size() / 2);

	for(int i = 0; i < (int)vertexInfo.size(); i++)
	{
		if(vertexInfo[i].checkIsBorder())
		{
			vertexInfo[i].flag = VF_BORDER;
			borderVerts.push_back(i);
		}
	}

	return borderVerts;
}

void Mesh::translate(Vec delta)
{
	this->translate(delta.x, delta.y, delta.z);
}

void Mesh::translate(double x, double y, double z)
{
	for( int i = 0; i < (int)vertex.size(); i++)
		vertex[i].add(x,y,z);

	vbo->setDirty(true);
}

void Mesh::translateVertices(const Vector<int> & vertices, Vec & delta)
{
	foreach(const int& i, vertices)
	{
		vertex[i].add(delta.x, delta.y, delta.z);
	}

	vbo->setDirty(true);
}

void Mesh::rotateVertices(const Vector<int> & vertices, const qglviewer::Quaternion & q, const Vec & pivot)
{
	foreach(const int& i, vertices)
	{
		vertex[i] = Point3D::RotateAround(vertex[i], pivot, q);
	}

	vbo->setDirty(true);
}

void Mesh::rotate(Vec axis, double angle)
{
	qglviewer::Quaternion q(axis, angle);

	for(int i = 0; i < (int)vertex.size(); i++)
		vertex[i].set(q.rotate(vertex[i]));

	vbo->setDirty(true);
}

void Mesh::scale(double factor)
{
	for(int i = 0; i < (int)vertex.size(); i++)
		vertex[i] *= factor;

	vbo->setDirty(true);
}

void Mesh::setDirtyVBO(bool state)
{
	vbo->setDirty(state);
}

void Mesh::rebuildVBO()
{
	vbo = new VBO(&this->vertex, &this->vNormal, &this->vColor, &this->face);

	this->isReady = true;
}

void Mesh::draw()
{
	if(this->isVisible && this->isReady)
	{
		// Enable alpha blending for transparent objects
		if(this->isTransparent)
		{
			glEnable(GL_BLEND); 
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		}

		// Change default shade to flat shading
		if(this->isFlatShade)
			glShadeModel(GL_FLAT);

		// use VBO (if available) or simple render method
		if(this->vbo && this->vbo->isVBOEnabled)
		{
			if(this->isDrawSmooth)		vbo->render_smooth();
			if(this->isDrawWireframe)	vbo->render_wireframe();
			if(this->isDrawVertices)	vbo->render_vertices();
			if(this->isDrawAsPoints)	this->drawSimplePoints(8);// vbo->render_as_points(); (ugly)
		}
		else
		{
			if(this->isDrawSmooth) this->drawSimple();
			if(this->isDrawWireframe) this->drawSimpleWireframe();
		}

		// Redo alpha render settings
		if(this->isTransparent)		glDisable(GL_BLEND);
		if(this->isFlatShade)		glShadeModel(GL_SMOOTH);

		if(this->isShowVertexNormals)
		{
			Vector<Vec> start, direction;
			for(int i = 0; i < (int)vertex.size(); i++){
				start.push_back(vertex[i]); direction.push_back(vNormal[i]);}
			SimpleDraw::DrawLineTick(start, direction);
		}

		if(this->isShowFaceNormals)
		{
			Vector<Vec> start, direction;
			for(int i = 0; i < (int)fNormal.size(); i++){
				start.push_back(f(i)->center()); direction.push_back(fNormal[i]);}
			SimpleDraw::DrawLineTick(start, direction, 0.5f, true, 1, 1, 0, 1);
		}
	}

	// Selected items
	double min_edge = DBL_MAX;

	foreach(int vi, selectedVertices.set)
	{
		// Compute it once
		if(min_edge == DBL_MAX)	{
			min_edge = minEdgeLengthAround(vi);
		}

		glEnable(GL_BLEND); 
		glColor4f(1,1,0,0.8f);
		SimpleDraw::DrawSphere(vertex[vi], min_edge * 0.20);
		glColor4f(1,1,1,1);
		glDisable(GL_BLEND);
	}

	// ===================
	// Debug: 
	for(StdMap<Point3D*, Color4>::iterator it = markers.begin(); it != markers.end(); it++)
	{
		Color4 * c = &it->second;
		SimpleDraw::IdentifyPoint(*it->first, c->r(), c->g(), c->b());
	} // draw markers

	for(Vector<Face>::iterator it = testFaces.begin(); it != testFaces.end(); it++)
		SimpleDraw::DrawTriangle(it->vec(0), it->vec(1), it->vec(2)); // test faces

	// test vertex
	Vector<Vec> testVecs;
	foreach(Point3D v, testVertex) testVecs.push_back(v);
	SimpleDraw::IdentifyPoints(testVecs); 

	for(Vector<Vector<Vec> >::iterator it = redPoints.begin(); it != redPoints.end(); it++)
		SimpleDraw::IdentifyConnectedPoints(*it, 1.0, 0.2f, 0.4f); // red points
	for(Vector<Vector<Vec> >::iterator it = bluePoints.begin(); it != bluePoints.end(); it++)
		SimpleDraw::IdentifyConnectedPoints(*it, 0.1f, 0.2f, 1.0); // blue points

	for(Vector<StdSet<Face*> >::iterator f = greenFaces.begin(); f != greenFaces.end(); f++){
		int i = 0; int N = f->size()/2;
		foreach(Face * face, *f){
			double alpha = (double)i / N;
			SimpleDraw::DrawTriangle(face->vec(0),face->vec(1),face->vec(2),0.2f,1,0.1f,alpha); i++;
		}
	} // green faces
	for(Vector<StdSet<Face*> >::iterator f = yellowFaces.begin(); f != yellowFaces.end(); f++){
		int i = 0; int N = f->size()/2;
		foreach(Face * face, *f){
			double alpha = (double)i / N;
			SimpleDraw::DrawTriangle(face->vec(0),face->vec(1),face->vec(2),1.0f,1.0f,0.1f,alpha); i++;
		}
	} // yellow faces

	SimpleDraw::IdentifyConnectedPoints(testPoints);// test points
	for(Vector<Line>::iterator it = testLines.begin(); it != testLines.end(); it++)	it->draw(); // draw test lines
	for(Vector<Plane>::iterator it = testPlanes.begin(); it != testPlanes.end(); it++)	it->draw(); // draw test planes

	if(selectedFace >= 0 && selectedFace < (int)face.size())
	{
		Face * face = faceIndexMap[selectedFace];
		glClear(GL_DEPTH_BUFFER_BIT);
		SimpleDraw::DrawTriangle(face->vec(0),face->vec(1),face->vec(2),1,1,1,1);

		SimpleDraw::IdentifyPoint(vertex[selectedVertex], 0,0,0,20);

		Umbrella u(&vertexInfo[selectedVertex]);
		Vector<Vec> points;

		for(int i = 0; i < (int)u.neighbor.size(); i++)
		{
			points.push_back(vec(u.neighbor[i]));
		}

		SimpleDraw::IdentifyConnectedPoints(points);
	}

	// === end debug
}

void Mesh::drawSimple()
{
	Vec f_normal, v1, v2, v3, n1, n2, n3;

	// This is useful for simple mesh constructions
	if((int)vNormal.size() != numberOfVertices())computeNormals();
	if((int)vColor.size() != numberOfVertices()) vColor = Vector<Color4>(vertex.size(), Color4());

	if(isFlatShade)
		glShadeModel(GL_FLAT);
	else
		glShadeModel(GL_SMOOTH);

	glBegin(GL_TRIANGLES);

	for(StdList<Face>::iterator f = face.begin(); f != face.end(); f++)
	{
		f_normal = fNormal[f->index];

		Color4 * color1 = &vColor[f->vIndex[0]];
		Color4 * color2 = &vColor[f->vIndex[1]];
		Color4 * color3 = &vColor[f->vIndex[2]];

		if(!isFlatShade){
			n1 = vNormal[f->vIndex[0]];
			n2 = vNormal[f->vIndex[1]];
			n3 = vNormal[f->vIndex[2]];
		}
		else
		{
			glNormal3dv(f_normal);
		}

		v1 = f->v[0]->vec();
		v2 = f->v[1]->vec();
		v3 = f->v[2]->vec();
		
		glColor4f(color1->r(), color1->g(), color1->b(), color1->a());
		if(!isFlatShade) glNormal3dv(n1);
		glVertex3dv(v1);

		glColor4f(color2->r(), color2->g(), color2->b(), color2->a());
		if(!isFlatShade) glNormal3dv(n2);
		glVertex3dv(v2);

		glColor4f(color3->r(), color3->g(), color3->b(), color3->a());
		if(!isFlatShade) glNormal3dv(n3);
		glVertex3dv(v3);
	}

	glEnd();
}

void Mesh::drawSimpleWireframe(bool useDefaultColor)
{
	glLineWidth(2.0); // default?

	glDisable(GL_LIGHTING);

	if(useDefaultColor)
		glColor3f(0, 0.75, 0);

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	Vec v1, v2, v3;

	glBegin(GL_TRIANGLES);
	for(StdList<Face>::iterator f = face.begin(); f != face.end(); f++)
	{
		v1 = f->v[0]->vec();
		v2 = f->v[1]->vec();
		v3 = f->v[2]->vec();

		glVertex3f(v1[0], v1[1], v1[2]);
		glVertex3f(v2[0], v2[1], v2[2]);
		glVertex3f(v3[0], v3[1], v3[2]);
	}
	glEnd();

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_LIGHTING);
}

void Mesh::drawSimplePoints(double pointSize)
{
	glPointSize(pointSize);
	glDisable(GL_BLEND);

	glBegin(GL_POINTS);
	for(int i = 0; i < (int)vertex.size(); i++)
	{
		if(vColor[i].a() > 0)
		{
			glColor4ubv(vColor[i]);
			glNormal3fv(vNormal[i]);
			glVertex3fv(vertex[i]);
		}
	}
	glEnd();
}

void Mesh::drawVertex(int i)
{
	glVertex3fv(vertex[i]);
}

void Mesh::drawEdge(Edge * e)
{
	glVertex3fv(vertex[e->vIndex[0]]);
	glVertex3fv(vertex[e->vIndex[1]]);
}

void Mesh::drawFace(Face * f)
{
	Vec v1 = f->v[0]->vec();
	Vec v2 = f->v[1]->vec();
	Vec v3 = f->v[2]->vec();
	Vec fn = fNormal[f->index];

	glNormal3fv(fn);
	glVertex3fv(v1);
	glVertex3fv(v2);
	glVertex3fv(v3);
}

void Mesh::drawUmbrella(Umbrella * u)
{
	glDisable(GL_LIGHTING);

	glPointSize(6.0);

	glBegin(GL_POINTS);
	glColor3f(1.0, 0, 0);
	for(int i = 0; i < (int)u->neighbor.size(); i++)
		drawVertex(u->neighbor[i]);
	glColor3f(0, 1.0, 0);
	drawVertex(u->index);
	glEnd();

	glBegin(GL_LINES);
	glLineWidth(3.0);
	glColor3f(0, 1.0, 1.0);
	foreach(HalfEdge h, u->halfEdge)
		drawEdge(&h.edge);
	glEnd();

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset( 0.5, 1.0 );

	/*glBegin(GL_TRIANGLES);
	glColor3f(0, 0, 1.0);
	foreach(Face * f, u->ifaces)
	drawFace(f);
	glEnd();*/

	glDisable(GL_POLYGON_OFFSET_FILL);

	glEnable(GL_LIGHTING);
}

void Mesh::drawVertexNames()
{
	for(int i = 0; i < (int)vertex.size(); i++)
	{
		glPushName(i);
		glBegin(GL_POINTS);
		glVertex3fv(vertex[i]);
		glEnd();
		glPopName();
	}
}

void Mesh::drawFaceNames()
{
	for(StdList<Face>::iterator f = face.begin(); f != face.end(); f++)
	{
		glPushName(f->index);
		glBegin(GL_TRIANGLES);
		glVertex3fv(f->vec(0));
		glVertex3fv(f->vec(1));
		glVertex3fv(f->vec(2));
		glEnd();
		glPopName();
	}
}

Mesh * Mesh::CloneSubMesh ( Vector<int> & facesIndex, bool isShallowClone, StdString newId)
{
	int numFaces = facesIndex.size();
	int numVertices = numFaces * 3;

	Mesh * clone = new Mesh();

	// String identification
	if(newId.size() > 0)
		clone->id = newId;
	else
		clone->id = this->id + "_clone";

	// Reserve known data structure sizes
	clone->vertex.reserve(numVertices);
	clone->vertexInfo.reserve(numVertices);

	int vCount = 0;
	int fCount = 0;
	int v1, v2, v3;

	// vIndexMap[ old index ] = new vertex index
	StdMap<int,int> vIndexMap;
	HashMap<int,bool> visisted;

	std::pair< std::map<int,int>::iterator , bool > r1, r2, r3;

	// Add selected faces to our sub-clone
	for(int i = 0; i < (int)facesIndex.size(); i++)
	{
		Face * face = this->faceIndexMap[facesIndex[i]];

		// For each vertex:
		// - Assign in map
		// - insert into mesh if didn't exists before
		r1 = vIndexMap.insert( PairInt(face->VIndex(0), vCount) );
		v1 = r1.first->second;
		if(r1.second){
			clone->addVertex( face->P(0)->vec(), v1 );
			vCount++;
		}
		r2 = vIndexMap.insert( PairInt(face->VIndex(1), vCount) );
		v2 = r2.first->second;
		if(r2.second){
			clone->addVertex( face->P(1)->vec(), v2 );
			vCount++;
		}
		r3 = vIndexMap.insert( PairInt(face->VIndex(2), vCount) );
		v3 = r3.first->second;
		if(r3.second){
			clone->addVertex( face->P(2)->vec(), v3 );
			vCount++;
		}

		// Add face
		clone->addFace(v1, v2, v3, fCount++);
	}

	if(!isShallowClone)
	{
		clone->clearColors();

		clone->computeNormals();
		clone->computeBounds();
		clone->createVBO();

		clone->isReady = true;
	}

	return clone;
}

void Mesh::addNoise(double scale)
{
	for( int i = 0; i < (int)vertex.size(); i++){
		double noiseX = (float)((rand() - (RAND_MAX / 2.0)) / RAND_MAX) * scale;
		double noiseY = (float)((rand() - (RAND_MAX / 2.0)) / RAND_MAX) * scale;
		double noiseZ = (float)((rand() - (RAND_MAX / 2.0)) / RAND_MAX) * scale;

		vertex[i].add(noiseX, noiseY, noiseZ);
	}

	this->computeNormals();
	this->computeBounds();
	
	if(octree)
		this->rebuildOctree();

	vbo->setDirty(true);
}

void Mesh::reassignFaces()
{
	for(StdList<Face>::iterator f = face.begin(); f != face.end(); f++)
	{
		f->v[0] = &this->vertex[f->vIndex[0]];
		f->v[1] = &this->vertex[f->vIndex[1]];
		f->v[2] = &this->vertex[f->vIndex[2]];
	}
}

void Mesh::refreshFaces(StdSet<Face *>& modifiedFaces)
{
	StdSet<VertexDetail*> modifiedNormals;

	int v1, v2, v3;
	Face * f;

	for(StdSet<Face*>::iterator it = modifiedFaces.begin(); it != modifiedFaces.end(); it++)
	{
		f = *it;

		v1 = f->vIndex[0];
		v2 = f->vIndex[1];
		v3 = f->vIndex[2];

		f->v[0] = &this->vertex[v1];
		f->v[1] = &this->vertex[v2];
		f->v[2] = &this->vertex[v3];

		// compute normals
		fNormal[f->index].set(f->normal());

		// also, later compute per vertex normals
		modifiedNormals.insert(&this->vertexInfo[v1]);
		modifiedNormals.insert(&this->vertexInfo[v2]);
		modifiedNormals.insert(&this->vertexInfo[v3]);
	}

	Normal v_normal;

	// Compute normals for modified vertices
	foreach(const VertexDetail * vd, modifiedNormals)
	{
		int i = vd->index;
		Vector<Face *> currFace = vertexInfo[i].ifaces;

		v_normal.setValue(0,0,0);

		for (int j = 0; j < (int)vertexInfo[i].ifaces.size(); j++)
			v_normal += fNormal[currFace[j]->index];

		vNormal[i] = Point3D(v_normal.unit());
	}
}

double Mesh::angleAroundVertex(int vIndex)
{
	double angle = 0;

	VertexDetail * vd = &vertexInfo[vIndex];
	Point3D * v = &vertex[vIndex];

	foreach(Face * f, vd->ifaces)
	{
		Edge edg = f->oppositeEdge(vIndex);

		angle += v->angleBetween(vertex[edg[0]], vertex[edg[1]]);
	}

	return angle;
}

int Mesh::vertexIndexClosest(const Vec& point)
{
	int closestIndex = -1;
	double minDist = DBL_MAX;

	for(int i = 0; i < (int)vertex.size(); i++)
	{
		double dist = (vertex[i] - point).norm();

		if(dist < minDist)
		{
			minDist = dist;
			closestIndex = i;
		}
	}

	return closestIndex;
}

void Mesh::clearAllVertexFlag()
{
#pragma omp parallel for
	for(int i = 0; i < (int)vertexInfo.size(); i++)
		vertexInfo[i].flag = VF_CLEAR;
}

void Mesh::clearAllFaceFlag()
{
	for(StdList<Face>::iterator f = face.begin(); f != face.end(); f++)
	{
		f->flag = FF_CLEAR;
	}
}

Umbrella Mesh::getUmbrella(int vertexIndex)
{
	return Umbrella(&vertexInfo[vertexIndex]);
}

void Mesh::getUmbrellas(Vector<Umbrella> & result)
{
	int N = vertexInfo.size();

	tempUmbrellas = Vector<Umbrella>(N);

	#pragma omp parallel for
	for(int i = 0; i < N; i++)
		tempUmbrellas[i] = Umbrella(&vertexInfo[i]);

	if(&tempUmbrellas != &result)
		result = tempUmbrellas;
}

void Mesh::getUmbrellas()
{
	getUmbrellas(tempUmbrellas);
}

HashMap<int, Vec> Mesh::getPoints()
{
	HashMap<int, Vec> result;

	for(int i = 0; i < (int)vertex.size(); i++)
		result[i] = vertex[i];

	return result;
}

int Mesh::getVertexIndexFromPos(const Vec& pos)
{
	int index = -1;

	for(int i = 0; i < (int)vertex.size(); i++)
	{
		if(pos.x == vertex[i].x && pos.y == vertex[i].y && pos.z == vertex[i].z)
		{
			index = i;
			break;
		}
	}

	return index;
}

StdSet<int> Mesh::getVerticesFromFaces(const Vector<int> & facesIndex)
{
	StdSet<int> verts;

	for(int i = 0; i < (int)facesIndex.size(); i++)
	{
		Face * face = this->faceIndexMap[facesIndex[i]];

		verts.insert(face->vIndex[0]);
		verts.insert(face->vIndex[1]);
		verts.insert(face->vIndex[2]);
	}

	return verts;
}

StdList<Face*> Mesh::getFacesFromIndices(const Vector<int> & facesIndex)
{
	StdList<Face*> result;

	for(int i = 0; i < (int)facesIndex.size(); i++)
		result.push_back(this->faceIndexMap[facesIndex[i]]);

	return result;
}

StdSet<Face*> Mesh::getFacesFromVertices(const StdSet<int> & vertexIndex, bool isExclusive)
{
	StdSet<Face*> result;

	for(StdSet<int>::const_iterator it = vertexIndex.begin(); it != vertexIndex.end(); it++)
	{
		VertexDetail * vd = &vertexInfo[*it];
		for(Vector<Face*>::iterator f = vd->ifaces.begin(); f != vd->ifaces.end(); f++)
		{
			if(!isExclusive)
			{
				result.insert(*f);
			}
			else
			{
				// Only include those faces with vertices in the input set
				if(SET_HAS(vertexIndex, (*f)->vIndex[0]) && 
					SET_HAS(vertexIndex, (*f)->vIndex[1]) && 
					SET_HAS(vertexIndex, (*f)->vIndex[2]))
				{
					result.insert(*f);
				}
			}
		}
	}

	return result;
}

Vector<Point3D> Mesh::getCopyPoints()
{
	return vertex;
}

const Vector<Point3D>& Mesh::getVertices() const
{
	return vertex;
}

const Vector<Normal>& Mesh::getNormals() const
{
	return vNormal;
}

void Mesh::setMeshPoints(const Vector<Point3D> & fromPoint)
{
	vertex = fromPoint;

	setDirtyVBO(true);
}

StdSet<int> Mesh::getConnectedPart(int vIndex)
{
	StdSet<int> partVerts;

	Vector<bool> visited(vertex.size(), false);
	Vector<int> * adj;

	Stack<int> stack;
	stack.push(vIndex);

	while(!stack.empty())
	{
		int cur = stack.top();
		stack.pop();

		adj = &tempUmbrellas[cur].neighbor;

		if(adj->size())
		{
			visited[cur] = true;
			partVerts.insert(cur);

			for(Vector<int>::iterator q = adj->begin(); q != adj->end(); q++) 
			{
				if (!visited[*q]) 
					stack.push(*q);
			}
		}
	}

	return partVerts;
}

Vector<int> Mesh::getFacesConnected( const Vector<int> & facesIndex, int firstFace)
{
	Vector<int> result;

	// Starting point
	if(firstFace < 0) firstFace = facesIndex.front();
	int firstVertex = faceIndexMap[firstFace]->vIndex[0];

	// Convert to vertices
	StdSet<int> vindices = getVerticesFromFaces(facesIndex);

	// Convert to connected graph
	Graph<int,float> g;
	foreach(int vi, vindices){
		foreach(int vj, vertexInfo[vi].adjacentVertices()){
			g.AddEdge(vi, vj, 1);
		}
	}

	// Explore starting from specified start
	StdSet<int> connectedVertices;
	g.explore(firstVertex, connectedVertices);

	// Add indices to a vector
	StdSet<Face*> connectedFaces = getFacesFromVertices(connectedVertices, false);
	foreach(Face * f, connectedFaces)
	{
		result.push_back(f->index);
	}

	return result;
}

void Mesh::removeAllFaces(const StdSet<int> & facesIndices)
{
	this->isReady = false;

	foreach(int fIndex, facesIndices)
	{
		faceIndexMap[fIndex]->unset();
	}

	StdSet<int> visited;
	HashMap<int,int> vIndexMap;

	int vindex = 0;

	Vector<Point3D> newVertices;
	Vector<Face> newFaces;

	for(StdList<Face>::iterator f = face.begin(); f != face.end(); f++)
	{
		if(f->VIndex(0) != -1 && f->VIndex(1) != -1 && f->VIndex(2) != -1)
		{
			for(int i = 0; i < 3; i++)
			{
				int vi = f->VIndex(i);

				if(!SET_HAS(visited, vi))
				{
					visited.insert(vi);
					vIndexMap[vi] = vindex++;
					newVertices.push_back(vertex[vi]);
				}
			}

			newFaces.push_back(*f);
		}
	}

	vertex.clear();
	vertexInfo.clear();

	int vIndex = 0;

	foreach(Point3D v, newVertices)
		addVertex(v, vIndex++);

	int fIndex = 0;

	face.clear();
	faceIndexMap.clear();

	foreach(Face f, newFaces)
		addFace(vIndexMap[f.vIndex[0]], vIndexMap[f.vIndex[1]], vIndexMap[f.vIndex[2]], fIndex++);
}

HoleStructure Mesh::getHoles()
{
	// Pre process: get border vertices and edges around vertices
	Vector<int> borderVertices = getBorderVertices();
	if(tempUmbrellas.size() < 1)	getUmbrellas(tempUmbrellas);

	// Hole structure
	HoleStructure hole;
	HashMap<int, int> belongsToHole;

	// fill belong map with -1
	for(int i = 0; i < (int)borderVertices.size(); i++)
		belongsToHole[borderVertices[i]] = -10;

	int currHole = numberOfVertices() - 1;
	int vertNotBelong = belongsToHole.size();

	while(vertNotBelong > 0)
	{
		// pick a vertex not belonging to a hole
		int v = -1;
		for(HashMap<int, int>::iterator it = belongsToHole.begin(); it != belongsToHole.end() ; it++)
		{
			if(it->second < 0)
			{
				v = it->first;
				break;
			}
		}
		if(v == -1)	break; // should not happen though

		// Declare a new hole
		currHole++;
		bool endOfLoop = false;

		// Add v to this hole
		belongsToHole[v] = currHole;
		hole[currHole].push_back(v);
		vertNotBelong--;

		// Find ring for the hole starting from v
		while( !endOfLoop )
		{
			PairInt n = tempUmbrellas[v].borderNeighbours();

			if(belongsToHole[n.first] == currHole && belongsToHole[n.second] == currHole)
			{
				endOfLoop = true;
				break;
			}

			if(belongsToHole[n.first] == currHole)
			{
				belongsToHole[n.second] = currHole;
				hole[currHole].push_back(n.second);

				v = n.second;
			}
			else
			{
				belongsToHole[n.first] = currHole;
				hole[currHole].push_back(n.first);

				v = n.first;
			}

			vertNotBelong--;
		}
	}

	return hole;
}

StdList<int> Mesh::getBoundry(int startIndex, bool rebuildUmbrellas)
{
	if(tempUmbrellas.size() < 1)
		getUmbrellas();

	Vector<Umbrella> * U = &tempUmbrellas;

	StdList<int> boundry;

	int start, next;

	start = startIndex;
	next = start;

	bool foundNext = false;

	if(rebuildUmbrellas)
		U->at(next) = Umbrella(&vertexInfo[next]);

	PairInt adj = U->at(next).borderNeighbours();

	int prev = adj.first;

	if(prev < 0)
		return boundry;

	boundry.push_back(prev);

	do
	{
		foundNext = false; // fail-safe

		boundry.push_back(next);

		if(rebuildUmbrellas)
			U->at(next) = Umbrella(&vertexInfo[next]);

		adj = U->at(next).borderNeighbours();

		if(adj.first == prev)
		{
			next = adj.second;	
			foundNext = true;
		}
		else if(adj.second == prev)	
		{
			next = adj.first;	
			foundNext = true;
		}

		prev = boundry.back();

	} while( foundNext && next != boundry.front() );

	return boundry;
}

StdSet<int> Mesh::visitFromBoundry(int boundryVertex, const StdSet<int>& border)
{
	StdSet<int> visited;
	std::queue<int> q;
	int vi	= boundryVertex;

	q.push(vi);
	visited.insert(vi);

	// First part
	while(!q.empty())
	{
		int i = q.front();
		q.pop();

		foreach(HalfEdge e, tempUmbrellas[i].halfEdge)
		{
			int j = e.edge.neighbor();

			// Check: not visited && not on inner border
			if(visited.find(j) == visited.end() && border.find(j) == border.end()) 
			{
				visited.insert(j);
				q.push(j);
			}
		}
	}

	return visited;
}

Face * Mesh::intersectRay( const Ray& ray, HitResult & hitRes )
{
	double minDist = DBL_MAX;
	Face *closestFace = NULL;

	double u = 0.0, v = 0.0;
	double actualMinDist = 0;

	for(StdList<Face>::iterator f = face.begin(); f != face.end(); f++)
	{
		f->intersectionTest(ray, hitRes, true);

		if(hitRes.hit && (abs(hitRes.distance) < minDist))
		{
			closestFace = &(*f);

			minDist = abs(hitRes.distance);
			actualMinDist = hitRes.distance;

			u = hitRes.u; 
			v = hitRes.v;
		}
	}

	if(closestFace)
	{
		// set 'hitRes' to that closest hit result
		hitRes.distance = actualMinDist;	
		hitRes.u = u;	hitRes.v = v;
		hitRes.hit = true;
	}
	else
	{
		hitRes.hit = false;
		hitRes.distance = DBL_MAX;
	}

	return closestFace;
}

StdList<BaseTriangle*> Mesh::facesListPointers()
{
	StdList<BaseTriangle*> result;

	for(StdList<Face>::iterator f = face.begin(); f != face.end(); f++)
	{
		result.push_back((BaseTriangle*)&(*f));
	}

	return result;
}

Face * Mesh::getBoundryFace( int vi1, int vi2 )
{
	return vertexInfo[vi1].sharedFace(vi2);
}

double Mesh::minEdgeLengthAround(int vi)
{
	double min_edge = DBL_MAX;

	StdSet<int> ne = vertexInfo[vi].adjacentVertices();

	foreach(int j, ne) min_edge = Min(min_edge, (vertex[vi]-vertex[j]).norm());

	return min_edge;
}

double Mesh::maxEdgeLengthAround(int vi)
{
	double max_edge = DBL_MIN;

	StdSet<int> ne = vertexInfo[vi].adjacentVertices();

	foreach(int j, ne) max_edge = Max(max_edge, (vertex[vi]-vertex[j]).norm());

	return max_edge;
}

StdSet<int> Mesh::getManifoldFaces( int startFace )
{
	StdSet<int> manifoldFaces;

	Stack<int> unvisitedVertices;
	IntSet visitedVertices;

	unvisitedVertices.push(faceIndexMap[startFace]->vIndex[0]);

	while(!unvisitedVertices.empty())
	{
		int vi = unvisitedVertices.top(); unvisitedVertices.pop();

		if(!vertexInfo[vi].isMissingBorder() && !visitedVertices.has(vi))
		{
			// add the adjacent faces
			for(int i = 0; i < (int)vertexInfo[vi].ifaces.size(); i++)
			{
				manifoldFaces.insert(vertexInfo[vi].ifaces[i]->index);
			}

			StdSet<int> adj = vertexInfo[vi].adjacentVertices();

			foreach(int j, adj)
			{
				if( !visitedVertices.has( j ) ) 
					unvisitedVertices.push(j);
			}

			// Mark as visisted
			visitedVertices.insert(vi);
		}
	}

	return manifoldFaces;
}

void Mesh::rebuildOctree()
{
	if(octree)
		delete octree;

	octree = new Octree();

	printf("Rebuilding Octree.."); CreateTimer(timer);

	octree->initBuild(facesListPointers(), 20);

	printf("Done (%d ms).\n", (int)timer.elapsed());
}

double Mesh::maxFaceArea()
{
	double maxArea = DBL_MIN;

	for(StdList<Face>::iterator f = face.begin(); f != face.end(); f++)
	{
		maxArea = Max(maxArea, f->area());
	}

	return maxArea;
}
