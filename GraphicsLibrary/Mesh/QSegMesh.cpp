#include "GUI/global.h"
#include "QSegMesh.h"
#include <fstream>
#include <set>
#include <map>
#include <QFile>
#include <QTextStream>
#include <QVector>
#include "Utility/SimpleDraw.h"

int turnOffSegments = 0;

QSegMesh::QSegMesh()
{
	isReady = false;
	isDrawAABB = false;
	upVec = Vec3d(0,0,1);
	radius = 1.0;
}

QSegMesh::QSegMesh( const QSegMesh& from )
{
	this->isReady = from.isReady;

	this->bbmin = from.bbmin;
	this->bbmax = from.bbmax;
	this->center = from.center;
	this->radius = from.radius;

	for (int i=0;i<(int)from.segment.size();i++)
	{
		QSurfaceMesh *seg_mesh = new QSurfaceMesh(*from.segment[i]);
		segment.push_back(seg_mesh);
	}

	this->upVec = from.upVec;
}

QSegMesh& QSegMesh::operator=( const QSegMesh& rhs )
{
	this->isReady = rhs.isReady;

	this->bbmin = rhs.bbmin;
	this->bbmax = rhs.bbmax;
	this->center = rhs.center;
	this->radius = rhs.radius;

	for (int i=0;i<(int)rhs.segment.size();i++)
	{
		QSurfaceMesh *seg_mesh = new QSurfaceMesh(*rhs.segment[i]);
		segment.push_back(seg_mesh);
	}

	this->upVec = rhs.upVec;

	return *this;
}

void QSegMesh::checkObjSegmentation ( QString fileName, QString segFilename)
{
	// Check if there is .seg file already
	std::fstream inF(qPrintable(segFilename), std::ios::in);
	if (inF) return;


	QFile file(fileName);
	file.open(QIODevice::ReadOnly | QIODevice::Text);
	QTextStream in(&file);

	QVector<int> faceGroups;
	bool collectingFaces = false;

	while (!in.atEnd()){
		QString line = in.readLine();

		// Possible face group is coming
		if(line.startsWith("g ")){
			collectingFaces = true;
			faceGroups.push_back(0);

			// Get name for this segment
			segmentName.push_back(line.replace("g ", "").trimmed());

			// If no name on file, give it one
			if(segmentName.back() == "")
				segmentName.back() = QString("Segment %1").arg(segmentName.size());
		}
		
		// Count faces
		if(line.startsWith("f ")){
			if(collectingFaces)
				faceGroups.back() += 1;
			else
				return;
		}
	}

	// Write segmentation file
	QFile segFile(segFilename);
	segFile.open(QIODevice::WriteOnly | QIODevice::Text);
	QTextStream out(&segFile);

	out << (int)faceGroups.size() << "\n";
	out << "labels ";
	for (int i=0;i<segmentName.size(); i++)
	{
		out<<segmentName[i]<<" ";
	}
	out << "\n";


	int offset = 0;

	for(int i = 0; i < (int)faceGroups.size(); i++){
		int fCount = faceGroups[i];
		if(i > 0) offset += faceGroups[i-1];

		for(int j = 0; j < fCount; j++)
			out << j + offset << " " << i << "\n";
	}
}

void QSegMesh::read( QString fileName )
{
	// Load entire mesh geometry
	QSurfaceMesh mesh;
	mesh.read(qPrintable(fileName));

	QString file_name = fileName;
	QString segFilename = file_name.replace(file_name.lastIndexOf('.') + 1, 3, "seg");

	// Check if file include segments "g Object"
	checkObjSegmentation(fileName, segFilename);

	// Load segmentation file
	std::ifstream inF(qPrintable(segFilename), std::ios::in);

	//if (mesh.n_faces() < 1 || !inF)
	if(turnOffSegments)
	{
		// Unsegmented mesh
		segment.push_back(new QSurfaceMesh(mesh));
	}
	else
	{
		int nbSeg;
		inF >> nbSeg;

		nbSeg = Max(1, nbSeg);

		segmentName.clear();
		std::string str;
		inF >> str;
		if (str == "labels")
		{
			for (int i=0;i<nbSeg;i++)
			{
				inF >> str;
				segmentName.push_back(str.c_str());
			}
		}
		else
		{// Fall back
			inF.seekg(-(int)str.size(), std::ios::cur);
		}

		std::vector<int> faceSeg(mesh.n_faces());
		int fid, sid;
		for (int i=0;i<(int)mesh.n_faces()&&inF;i++)
		{
			inF >> fid >> sid;	
			faceSeg[fid] = sid;
		}
		inF.close();

		// Create segments
		for (int i=0;i<nbSeg;i++)
		{
			segment.push_back(new QSurfaceMesh());
		}


		// Create unique vertex set for each segment
        std::vector< std::set <Surface_mesh::Vertex> > segVertices(nbSeg);
		Surface_mesh::Face_iterator fit, fend = mesh.faces_end();
		Surface_mesh::Vertex_around_face_circulator fvit;	

		for (fit = mesh.faces_begin(); fit!=fend; ++fit)
		{
			Surface_mesh::Face f = fit;
			int sid = faceSeg[f.idx()];

			fvit = mesh.vertices(fit);	
			segVertices[sid].insert(fvit);
			segVertices[sid].insert(++fvit);
			segVertices[sid].insert(++fvit);
		}

		// Add Vertices to each segment	
        std::vector< std::map <Surface_mesh::Vertex, Surface_mesh::Vertex> > segVerMap(nbSeg);
		Surface_mesh::Vertex_property<Point>  points   = mesh.vertex_property<Point>("v:point");

		for (int i=0;i<nbSeg;i++)
		{
			std::set<Surface_mesh::Vertex>::iterator vit, vend = segVertices[i].end();
			int j = 0;
			for (vit=segVertices[i].begin(); vit!=vend; vit++, j++)
			{
				segment[i]->add_vertex(mesh.getVertexPos(*vit));

				segVerMap[i].insert(std::make_pair(*vit, mesh.getVertex(j)));  // Create a new index for each vertex
			}
		}

		// Add Faces to each segment
		std::vector<Surface_mesh::Vertex>  vertices(3);
		for (fit = mesh.faces_begin(); fit!=fend; ++fit)
		{
			Surface_mesh::Face f = fit;
			int sid = faceSeg[f.idx()];

			fvit = mesh.vertices(fit);
			vertices[0] = segVerMap[sid][fvit];
			vertices[1] = segVerMap[sid][++fvit];
			vertices[2] = segVerMap[sid][++fvit];

			segment[sid]->add_face(vertices);
		}

	}

	// Clear the empty segment
	for (std::vector<QSurfaceMesh*>::iterator itr=segment.begin(); itr!=segment.end(); )
	{
		if (!(*itr)->n_vertices())
		{
			itr = segment.erase(itr);
		}
		else
			itr++;
	}

	// Build up
	build_up();
}

void QSegMesh::saveObj( QString fileName )
{
	FILE * outF = fopen (qPrintable(fileName) , "w");

	uint faceOffset = 1;

	for (uint i = 0; i < nbSegments(); i++)
	{
		std::vector<Point> segPoints = segment[i]->clonePoints();
		std::vector<uint> segFaces = segment[i]->cloneTriangleIndices();

		foreach(Point p, segPoints)
			fprintf(outF, "v %f %f %f\n", p.x(), p.y(), p.z());
		fprintf(outF,"# %u vertices\n\n\n", segPoints.size()); // info

		fprintf(outF,"g %s\n", qPrintable(segment[i]->objectName()));

		for(int vi = 0; vi < (int)segFaces.size(); vi += 3)
		{
			uint v0 = segFaces[vi + 0], v1 = segFaces[vi + 1], v2 = segFaces[vi + 2];

			fprintf(outF, "f %u %u %u\n", faceOffset + v0, faceOffset + v1, faceOffset + v2);
		}
		fprintf(outF,"# %u faces\n\n\n", segment[i]->n_faces()); // info

		faceOffset += segment[i]->n_vertices();
	}

	fclose(outF);
}

void QSegMesh::insertCopyMesh(QSurfaceMesh * newSegment)
{
	this->segment.push_back(new QSurfaceMesh(*newSegment));
}

void QSegMesh::build_up()
{
	computeBoundingBox();
	moveCenterToOrigin();
	normalize();
	computeBoundingBox();

	update_face_normals();
	update_vertex_normals();

	for (int i = 0; i < (int)nbSegments();i++)
	{
		segment[i]->buildUp();
	}
	
	setColorVertices();

	printf("Segments loaded: %d \n", nbSegments());

	isReady = true;
}

void QSegMesh::update_face_normals()
{
	for (int i = 0; i < (int)nbSegments();i++)
		segment[i]->update_face_normals();
}

void QSegMesh::update_vertex_normals()
{
	for (int i = 0; i < (int)nbSegments();i++)
		segment[i]->update_vertex_normals();
}


void QSegMesh::computeBoundingBox()
{
	// compute bounding box
	bbmin = Point( FLT_MAX,  FLT_MAX,  FLT_MAX);
	bbmax = Point(-FLT_MAX, -FLT_MAX, -FLT_MAX);	
	
	for (int i = 0; i < (int)nbSegments();i++)
	{
		if (!segment[i]->isVisible) continue;

		Surface_mesh::Vertex_property<Point> points = segment[i]->vertex_property<Point>("v:point");
		Surface_mesh::Vertex_iterator vit, vend = segment[i]->vertices_end();

		for (vit = segment[i]->vertices_begin(); vit != vend; ++vit)
		{
			bbmin.minimize(points[vit]);
			bbmax.maximize(points[vit]);
		}		
	}
	
	center = (bbmin + bbmax) * 0.5f;
	radius = (bbmax - bbmin).norm() * 0.5f;
}

void QSegMesh::moveCenterToOrigin()
{
	if (!MOVE_CENTER_TO_ORIGIN) return;

	for (int i = 0; i < (int)nbSegments();i++)
	{
		Surface_mesh::Vertex_property<Point> points = segment[i]->vertex_property<Point>("v:point");
		Surface_mesh::Vertex_iterator vit, vend = segment[i]->vertices_end();

		for (vit = segment[i]->vertices_begin(); vit != vend; ++vit)
		{
			if (!segment[i]->is_deleted(vit))
				points[vit] -= center;
		}
	}
}

void QSegMesh::setColorVertices( double r, double g, double b, double a )
{
	bool assignRandomColors = false;
	std::vector<Vec4d> randomColors;

	if(assignRandomColors)
		randomColors = SimpleDraw::RandomColors(segment.size());
	else
	{
		randomColors = std::vector<Vec4d>(nbSegments(), Vec4d(r,g,b,a));
	}

	for (int i=0;i<(int)segment.size();i++)
	{
		segment[i]->setColorVertices(randomColors[i]);
	}
}

uint QSegMesh::nbSegments()
{
	return segment.size();
}

QSurfaceMesh* QSegMesh::operator[]( uint i )
{
	return segment[i];
}

QSurfaceMesh* QSegMesh::getSegment( uint i )
{
	return segment[i];
}

QSurfaceMesh* QSegMesh::getSegment( QString sid )
{
	for(int i = 0; i < segmentName.size(); i++)
		if(segmentName[i] == sid)
			return segment[i];
	return NULL;
}

std::vector<QSurfaceMesh*> QSegMesh::getSegments()
{
	return segment;
}

void QSegMesh::simpleDraw( bool isColored /*= true*/, bool isDots /*= false*/ )
{
	// Render mesh regularly (inefficient)
	for (int i = 0;i < (int)segment.size(); i++)
	{
		if (segment[i]->isVisible)
			segment[i]->simpleDraw(isColored, isDots);
	}
}

void QSegMesh::drawFacesUnique()
{
	uint offset = 0;
	for (int i=0;i<(int)segment.size();i++)
	{
		if (segment[i]->isVisible)
		{
			segment[i]->drawFacesUnique(offset);
			offset += segment[i]->n_faces();
		}
	}
}

void QSegMesh::drawDebug()
{
	for (int i=0;i<(int)segment.size();i++)
		segment[i]->drawDebug();

	if(isDrawAABB)
		drawAABB();
}

void QSegMesh::setObjectName( const QString &name )
{
	QObject::setObjectName(name);

	// For single objects, broken 'seg' files
	int i = 0;
	while(segmentName.size() < (int)segment.size())
		segmentName.push_back(name + QString("-%1").arg(i++));

	for (int i=0;i<(int)segment.size();i++)
	{
		segment[i]->setObjectName(segmentName[i]);
	}
}

uint QSegMesh::nbVertices()
{
	uint num = 0;

	for (int i=0;i<(int)segment.size();i++)
		num += segment[i]->n_vertices();

	return num;
}


uint QSegMesh::nbFaces()
{
	uint num = 0;

	for (int i=0;i<(int)segment.size();i++)
		num += segment[i]->n_faces();

	return num;
}

std::vector<uint> QSegMesh::vertexIndicesAroundFace( uint fid )
{
	uint sid;
	uint fid_local;
	global2local_fid(fid, sid, fid_local);

	std::vector<uint> vertices = segment[sid]->vertexIndicesAroundFace(fid_local);
	uint offset = fid - fid_local;

	for (int i=0;i<(int)vertices.size();i++)
		vertices[i] += offset;

	return vertices;
}

void QSegMesh::global2local_fid( uint fid, uint& sid, uint& fid_local )
{
	if (fid >= nbFaces()) return;

	uint offset = 0;
	int i=0;
	for (;i<(int)segment.size();i++)
	{


		offset += segment[i]->n_faces();

		if (fid < offset)
		{
			offset -= segment[i]->n_faces();
			break;
		}
	}

	sid = i;
	fid_local = fid - offset;
}

Point QSegMesh::getVertexPos( uint vid )
{
	uint sid;
	uint vid_local;
	global2local_vid(vid, sid, vid_local);

	return segment[sid]->getVertexPos(vid_local);
}

void QSegMesh::global2local_vid( uint vid, uint& sid, uint& vid_local )
{
	uint offset = 0;
	int i=0;
	for (;i<(int)segment.size();i++)
	{
		if (!segment[i]->isVisible)
			continue;

		offset += segment[i]->n_vertices();

		if (vid < offset)
		{
			offset -= segment[i]->n_vertices();
			break;
		}
	}

	sid = i;
	vid_local = vid - offset;
}

void QSegMesh::setVertexColor( uint vid, const Color& newColor )
{
	uint sid;
	uint vid_local;
	global2local_vid(vid, sid, vid_local);

	segment[sid]->setVertexColor(vid_local, newColor);
}

void QSegMesh::normalize()
{
	if (!NORMALIZE_MESH) return;

	for (uint i = 0; i < nbSegments();i++)
	{
		Surface_mesh::Vertex_property<Point> points = segment[i]->vertex_property<Point>("v:point");
		Surface_mesh::Vertex_iterator vit, vend = segment[i]->vertices_end();

		for (vit = segment[i]->vertices_begin(); vit != vend; ++vit)
		{
			if (!segment[i]->is_deleted(vit))
				points[vit] = points[vit] / radius;
		}
	}

	scaleFactor = radius;
}

void QSegMesh::rotateAroundUp( double theta )
{
	for (uint i = 0; i < nbSegments();i++)
		segment[i]->rotateAroundUp(theta);
}

void QSegMesh::rotateUp( Vec3d to )
{
	double theta = acos(dot(upVec, to));
	if(theta == 0)	return;
	Vec3d axis = cross(upVec, to).normalized();

	for (uint i = 0; i < nbSegments();i++)
	{
		Vec3d oldCenter = segment[i]->center;
		Vec3d newCenter = ROTATED_VEC(oldCenter, theta, axis);
		Vec3d delta = newCenter - oldCenter;
		segment[i]->translate(delta);

		segment[i]->rotateUp(to);
	}

	upVec = to;
}

uint QSegMesh::segmentIdOfVertex( uint vid )
{
	uint sid;
	uint vid_local;
	global2local_vid(vid, sid, vid_local);

	return sid;
}


//		  4-----------7                     Z
//		 /|          /|                     ^   /
//		5-+---------6 |                     |  / 
//		| |         | |                     | /
//		| |         | |                     |/     
//		| 0---------+-3            ---------+-------> Y 
//		|/          |/                     /|
//		1-----------2                     / |
//								        ↙  |
//	                                    X
void QSegMesh::drawAABB()
{
	Vec3d X(0, 0, 0), Y(0, 0, 0), Z(0, 0, 0);
	X[0] = bbmax[0]-bbmin[0];
	Y[1] = bbmax[1]-bbmin[1];
	Z[2] = bbmax[2]-bbmin[2];

	std::vector< Vec3d > P(8);
	P[0] = bbmin;
	P[1] = P[0] + X;
	P[2] = P[1] + Y;
	P[3] = P[0] + Y;

	P[4] = P[0] + Z;
	P[5] = P[1] + Z;
	P[6] = P[2] + Z;
	P[7] = P[3] + Z;

	SimpleDraw::IdentifyLine(P[0], P[1], Color(0, 1, 1,1), false);
	SimpleDraw::IdentifyLine(P[1], P[2], Color(0, 1, 1,1), false);
	SimpleDraw::IdentifyLine(P[2], P[3], Color(0, 1, 1,1), false);
	SimpleDraw::IdentifyLine(P[3], P[0], Color(0, 1, 1,1), false);

	SimpleDraw::IdentifyLine(P[4], P[5], Color(0, 1, 1,1), false);
	SimpleDraw::IdentifyLine(P[5], P[6], Color(0, 1, 1,1), false);
	SimpleDraw::IdentifyLine(P[6], P[7], Color(0, 1, 1,1), false);
	SimpleDraw::IdentifyLine(P[7], P[4], Color(0, 1, 1,1), false);

	SimpleDraw::IdentifyLine(P[0], P[4], Color(0, 1, 1,1), false);
	SimpleDraw::IdentifyLine(P[1], P[5], Color(0, 1, 1,1), false);
	SimpleDraw::IdentifyLine(P[2], P[6], Color(0, 1, 1,1), false);
	SimpleDraw::IdentifyLine(P[3], P[7], Color(0, 1, 1,1), false);
}

QSurfaceMesh * QSegMesh::flattenMesh()
{
	QSurfaceMesh * flattend = new QSurfaceMesh;
	Surface_mesh::Vertex v0, v1, v2;

	flattend->id = objectName();
	flattend->setObjectName(objectName());

	// Add geometry of all segments together
	int offset = 0;

	for (uint i = 0; i < nbSegments();i++)
	{
		Surface_mesh::Vertex_property<Point> points = segment[i]->vertex_property<Point>("v:point");
		Surface_mesh::Vertex_iterator vit, vend = segment[i]->vertices_end();

		// Add vertices
		for (vit = segment[i]->vertices_begin(); vit != vend; ++vit)
			flattend->add_vertex(points[vit]);

		Surface_mesh::Face_iterator fit, fend = segment[i]->faces_end();

		// Add faces
		for (fit = segment[i]->faces_begin(); fit!=fend; ++fit)
		{
			std::vector<uint> v = segment[i]->vertexIndicesAroundFace(Surface_mesh::Face(fit).idx());
			
			v[0] += offset;
			v[1] += offset;
			v[2] += offset;

			flattend->add_triangle(Surface_mesh::Vertex(v[0]), Surface_mesh::Vertex(v[1]), Surface_mesh::Vertex(v[2]));
		}

		// offset for next segment
		offset += segment[i]->n_vertices();
	}

	return flattend;
}

