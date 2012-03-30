#pragma once
#include <QObject>
#include <QString>
#include "Utility/Macros.h"

#include "GraphicsLibrary/Mesh/SurfaceMesh/Surface_mesh.h"

class QSurfaceMesh : public QObject, public Surface_mesh
{
	Q_OBJECT

public:
	QSurfaceMesh();
	QSurfaceMesh(const QSurfaceMesh& from);
	QSurfaceMesh& operator=(const QSurfaceMesh& rhs);

	std::vector<Vertex_iterator> vertex_array;
	std::vector<Face_iterator> face_array;

	void assignVertexArray();
	void assignFaceArray();

	std::vector<uint> vertexIndicesAroundFace( uint f_id );	
	Point getVertexPos( uint v_id );
	Point getVertexPos( const Vertex v );
	void setVertexPos( const Vertex v, Point newPos);
	Vertex getVertex( uint v_id);
	Face getFace( uint f_id);

	std::set<uint> faceIndicesAroundVertex(const Vertex& v);
	std::set<uint> vertexIndicesAroundVertex(const Vertex& v);

	void computeBoundingBox();
	void moveCenterToOrigin();
	void translate(Vec3d delta);

	void rotateUp( Vec3d to);
	Vec3d upVec;
	void rotateAroundUp(double theta);

	double getAverageEdgeLength();
	double averageEdgeLength;

	double volume();
	double normalize();
	double scalingFactor;

	Point closestPointVertices(const Point & p);
	double closestDistancePointVertices(const Point & p);
	uint closestVertex( const Point & p );
	Point closestPointFace(Face f, const Point & p);

	std::vector<Point> clonePoints();
	void setFromPoints(const std::vector<Point>& fromPoints);
	void setFromNormals( const std::vector<Normal>& fromNormals );
	std::vector< std::pair<Point, Point> > cloneEdges();

	std::vector<Normal> cloneFaceNormals();

	void drawFaceNames();
	void drawFacesUnique(uint offset);
	void drawDebug();
	void simpleDraw(bool isColored = true, bool isDots = false);
	void simpleDrawWireframe();

	void setColorVertices(double r = 1.0, double g = 1.0, double b = 1.0, double a = 1.0);
	void setColorVertices( Vec4d color );
	void setVertexColor( uint v_id, const Color& newColor );

	void collectEnoughRings(Vertex v, const size_t min_nb, std::vector <Vertex>& all);
	void resetVistedVertices(std::vector <Vertex>& all);
	void resetVistedVertices(uint toState = false); // for entire mesh

	void addNoise(double delta);
	void push( Vec3d from, Vec3d to, double falloff );

	// Load the mesh from file
	void read( const std::string& filename, bool isBuildUp = true);
	void writeObj(const std::string& filename);

	// Build up
	void buildUp();

	// Properties
	bool isReady;
	Point bbmin, bbmax, center;
	Scalar radius;

	// Debug items
	std::vector<Point> debug_points, debug_points2, debug_points3;
	std::vector< std::vector<Point> > debug_lines, debug_lines2, debug_lines3;
	bool isDrawBB;

	// Face Utility
	Vec3d fn( Face f );
	Vec3d faceCenter ( Face f );
	double faceArea( Face f );
	std::vector<Vec3d> facePoints( Face f );
	std::vector<uint> faceVerts( const Face& f );
	Vec3d getBaryFace( Face f, double U, double V );
	void fillTrianglesList();
	std::vector<unsigned int> cloneTriangleIndices();

	std::vector<unsigned int> triangles, edges;

	QString id;
	bool isVisible;

	std::vector<Vec3d> specialPnts;

private:
	bool isDirty;
};
