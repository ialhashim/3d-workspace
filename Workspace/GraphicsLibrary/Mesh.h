#pragma once

// Basic Structures
#include "Point.h"
#include "Face.h"
#include "VertexDetail.h"
#include "Edge.h"
#include "HalfEdge.h"
#include "Umbrella.h"

#include "Color4.h"
#include "Line.h"
#include "Plane.h"
#include "VBO.h"
#include "Octree.h"
#include "HashTable.h"

typedef std::map<int, Vector<int> > HoleStructure;

class Mesh
{
protected:

	// Vertices
	Vector<Point3D> vertex;
	Vector<VertexDetail> vertexInfo;

	// Faces
	StdList<Face> face;
	HashMap<int, Face* > faceIndexMap;

	// Normal & Color
	Vector<Normal> vNormal;
	Vector<Normal> fNormal;
	Vector<Color4> vColor;

	// Vertex Buffer Object
	VBO * vbo;

public:
	Mesh(int expectedNumVerts = 0);
	~Mesh();
	Mesh(const Mesh& fromMesh);					// Copy constructor
	Mesh& operator= (const Mesh& fromMesh);		// Assignment operator

	void mergeWith(const Mesh& other);

	// LOAD/SAVE MESH DATA FROM FILE
	void loadFromFile(const char* fileName);
	void saveToFile(const char* fileName = NULL);

	void loadFromFileOFF(const char* fileName);
	void loadFromFileOBJ(const char* fileName);

	bool isReady;

	// INTERSECTIONS
	Octree *octree;
	void rebuildOctree();

	StdString id;

	// Radius, Center, and Bounds
	double radius;
	Point3D center;
	Point3D minBound;
	Point3D maxBound;
	double normalize_scale;
	double maxScale;

	// Vertex & Face count
	int numberOfVertices();
	int numberOfFaces();

	// NORMALS & BOUNDS
	void computeNormals();
	void computeBounds();
	void moveToCenter();
	void normalizeScale();
	double computeVolume();
	double computeRadius();
	Vec computeCenter();

	Vec computeVNormalAt(int vi);

	// SIMPLE ELEMENTS CREATION
	void addVertex(double x, double y, double z, int index);
	void addVertex(Vec v, int index);
	void addFace(int v1, int v2, int v3, int index, bool forceOrientation = false);

	// ACCESSORS
	inline Point3D * v(int index)		{return &vertex[index];}		// Point3D pointer
	inline Normal * n(int index)		{return &vNormal[index];}		//  -Normal
	inline Color4 * vc(int index)		{return &vColor[index];}		//  -Color
	inline VertexDetail * vd(int index) {return &vertexInfo[index];}	//  -Detail
	inline Vec vec(int index)			{return vertex[index];}			//  -Position
	inline Point3D& ver(int index)		{return vertex[index];}			//  -Position (by reference)
	inline Face * f(int index)			{return faceIndexMap[index];}	// Face pointer
	inline Normal * fn(int index)		{return &fNormal[index];}		//  -Normal
	inline Umbrella * u(int index)		{return &tempUmbrellas[index];}	// Umbrella pointer

	bool withID(StdString ID);

	// Umbrellas
	Umbrella getUmbrella(int Point3DIndex);
	void getUmbrellas(Vector<Umbrella> & result);
	void getUmbrellas();
	Vector<Umbrella> tempUmbrellas;

	double angleAroundVertex(int vIndex);
	int vertexIndexClosest(const Vec& point);

	// ACCESS OPERATIONS
	StdList<Face> * facesList() { return &face; }
	StdList<BaseTriangle*> facesListPointers();
	void removeAllFaces(const StdSet<int> & facesIndices);
	HashMap<int, Vec> getPoints();

	// Access to vertices and normals
	const Vector<Point3D>& getVertices() const;
	const Vector<Normal>& getNormals() const;

	StdSet<int> getVerticesFromFaces(const Vector<int> & facesIndex);
	StdList<Face*> getFacesFromIndices(const Vector<int> & facesIndex);
	StdSet<Face*> getFacesFromVertices(const StdSet<int> & Point3DIndex, bool isExclusive = false);
	Vector<int> getFacesConnected(const Vector<int> & facesIndex, int firstFace = -1);
	Vector<Point3D> getCopyPoints();
	void setMeshPoints(const Vector<Point3D> & fromPoint);

	// VERTEX-BASED OPERATIONS
	Face * getBoundryFace(int vi1, int vi2);
	double minEdgeLengthAround(int vi);
	double maxEdgeLengthAround(int vi);
	int getVertexIndexFromPos(const Vec& pos);

	// FACE-BASED OPERATIONS
	double maxFaceArea();

	// HOLE OPERATIONS
	StdSet<int> getConnectedPart(int vIndex);
	HoleStructure getHoles();
	StdList<int> getBoundry(int vIndex, bool rebuildUmbrellas = false);
	StdSet<int> visitFromBoundry(int boundryPoint3D, const StdSet<int>& border);
	StdSet<int> getManifoldFaces(int startFace);

	// MODIFIERS
	void scale(double factor);
	void rotate(Vec axis, double angle);
	void translate(Vec delta);
	void translate(double x, double y, double z);
	void translateVertices(const Vector<int> & vertices, Vec & delta);
	void rotateVertices(const Vector<int> & vertices, const qglviewer::Quaternion & q, const Vec & pivot = Vec());
	void setColor(const Color4 & color);
	void setColor(int r, int g, int b, int a = 255);
	void setColor(Vector<int> & vertices, int r, int g, int b, int a = 255);
	void setColorFace(int faceIndex, int r, int g, int b, int a = 255);
	void setColorFaces(Vector<int> & faces, int r, int g, int b, int a = 255);
	void colorSplit(Vec & pos, Vec & dir);
	void clearColors();
	void addNoise(double scale);

	void reassignFaces();
	void refreshFaces(StdSet<Face *>& modifiedFaces);

	Face * intersectRay(const Ray& ray, HitResult & hitRes);

	// Vertex FLAGS
	void clearAllVertexFlag();
	void flagBorderVertices();
	Vector<int> getBorderVertices();
	void colorMarkBorders(int r, int g, int b, int a = 255);

	// FACE FLAGS
	void clearAllFaceFlag();

	// VBO OPERATIONS
	void createVBO();
	void updateVBO();
	void setDirtyVBO(bool state = true);
	void rebuildVBO();

	// RENDERING
	void draw();
	void drawSimple(); // if VBO not supported
	void drawSimpleWireframe(bool useDefaultColor = true);
	void drawSimplePoints(double pointSize = 6.0f);
	void drawVertex(int i);
	void drawEdge(Edge * e);
	void drawFace(Face * f);
	void drawUmbrella(Umbrella * u);

	// SELECTION
	void drawVertexNames();
	void drawFaceNames();

	IntSet selectedVertices;
	IntSet selectedFaces;

	// RENDERING FLAGS
	bool isVisible;
	bool isDrawSmooth;
	bool isDrawWireframe;
	bool isDrawVertices;
	bool isTransparent;
	bool isFlatShade;
	bool isDrawAsPoints;
	bool isShowVertexNormals;
	bool isShowFaceNormals;

	// DEBUG PROPERTIES
	StdMap<Point3D*, Color4> markers;
	Vector<Face> testFaces;
	Vector<Point3D> testVertex;
	Vector<Vec> testPoints;
	Vector<Line> testLines;
	Vector<Plane> testPlanes;
	Vector<Vector<Vec> > redPoints;
	Vector<Vector<Vec> > bluePoints;
	Vector<StdSet<Face*> > greenFaces;
	Vector<StdSet<Face*> > yellowFaces;
	int selectedVertex, selectedFace;

	// OTHER CONTAINERS FOR MISC DATA
	StdMap<StdString, double> item_flt;
	StdMap<StdString, int> item_int;
	StdMap<StdString, bool> item_bool;
	StdString opendFileName;

	// SUB MESH
	Mesh * CloneSubMesh(Vector<int> & facesIndex, bool isShallowClone = false, StdString newId = "");

	// My friends
	friend class Smoother;
	friend class Slicer;
	friend class Curvature;
	friend class Sampler;
	friend class Voxeler;
};
