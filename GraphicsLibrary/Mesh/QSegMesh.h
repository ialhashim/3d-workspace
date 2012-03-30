#pragma once
#include <QObject>
#include <QString>
#include <QVector>
#include <QMap>
#include "QSurfaceMesh.h"
#include <vector>

class QSegMesh : public QObject
{
	Q_OBJECT

public:
	QSegMesh();
	QSegMesh(const QSegMesh& from);
	QSegMesh& operator=(const QSegMesh& rhs);

	// Face, vertex
	uint nbVertices();
	uint nbFaces();
	std::vector<uint> vertexIndicesAroundFace(uint fid);
	Point getVertexPos( uint vid );
	void setVertexColor( uint vid, const Color& newColor );
	void global2local_fid( uint fid, uint& sid, uint& fid_local );
	void global2local_vid(uint vid, uint& sid, uint& vid_local);

	// Get segment
	QSurfaceMesh* operator [] (uint i);
	QSurfaceMesh* getSegment(uint i);
	QSurfaceMesh* getSegment( QString sid );
	std::vector<QSurfaceMesh*> getSegments();
	uint nbSegments();
	uint segmentIdOfVertex( uint vid );

	// Draw
	void simpleDraw(bool isColored = true, bool isDots = false);
	void drawFacesUnique();
	void drawDebug();
	void drawAABB();

	// Load the mesh from file
	void read(QString fileName);

	// Save the mesh
	void saveObj(QString fileName);

	// Build up the mesh
	void build_up();
	void moveCenterToOrigin();
	void computeBoundingBox();
	void setColorVertices( double r = 1.0, double g = 1.0, double b = 1.0, double a = 1.0);
	void update_face_normals();
	void update_vertex_normals();
	void normalize();
	void insertCopyMesh(QSurfaceMesh * newSegment);
	QSurfaceMesh * flattenMesh();

	// Basic rotations
	void rotateUp(Vec3d to);
	Vec3d upVec;
	void rotateAroundUp(double theta);

	// Properties
	bool isReady;
	bool isDrawAABB;
	Point bbmin, bbmax, center;
	Scalar radius;
	Scalar scaleFactor;

	// Set global unique name for this and all its segments
	void setObjectName(const QString &name);
	QVector<QString> segmentName;

	// Other auxiliary data
        QMap<QString, std::vector<double> > data1D;
        QMap<QString, std::vector< std::vector<double> > > data2D;

private:
	std::vector<QSurfaceMesh*> segment;
	
	// This is useful for segmented OBJs
	void checkObjSegmentation ( QString fileName, QString segFilename);

};
