#pragma once

#include "Point.h"
#include "Face.h"
#include "Color4.h"

class VBO
{
	unsigned int vertex_vbo_id;
	unsigned int normal_vbo_id;
	unsigned int color_vbo_id;
	unsigned int faces_id;

public:
	Vector<Point3D> * vertices;
	Vector<Normal> * normals;
	Vector<Color4> * colors;
	Vector<Index> * indices;

	VBO(Vector<Point3D> * v, Vector<Normal> * n, Vector<Color4> * c, StdList<Face> * f);
	~VBO();

	void update();
	void update_vbo(Index *vbo, int vbo_size, const void *vbo_data);
	void update_ebo(Index *ebo, int ebo_size, const void *ebo_data);
	void free_vbo(Index vbo);

	// Rendering Vertex Buffer Object (VBO)
	void render_smooth(bool dynamic = false);
	void render_wireframe(bool dynamic = false);
	void render_vertices(bool dynamic = false);
	void render_as_points(bool dynamic = false);

	bool isDirty;
	inline void setDirty(bool state) {isDirty = state;}

	bool isVBOEnabled;
};

