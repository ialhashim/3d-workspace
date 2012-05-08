#pragma once

#include "Utility/Macros.h"

// Surface-mesh
#include "Vector.h"

typedef Vec3d PointType;
typedef Vec3d NormalType;
typedef Vec4d ColorType;

enum RENDER_MODE{ RENDER_WIREFRAME, RENDER_POINT, RENDER_REGULAR};

class VBO
{
	unsigned int vertex_vbo_id;
	unsigned int Normalvbo_id;
	unsigned int color_vbo_id;
	unsigned int faces_id;

public:
	std::string objectId;
	int vCount;

	const PointType * vertices;
	const NormalType * normals;
	const ColorType * colors;
	StdVector<uint> indices;

	VBO(){ isVBOEnabled = false; isDirty = false; isReady = false; };
	VBO( unsigned int vert_count, const PointType * v, const NormalType * n, const ColorType * c, StdVector<uint> faces );

	void free_vbo(uint vbo);
	~VBO();

	static bool isVBOSupported();

	void update();
	void update_vbo(uint *vbo, int vbo_size, const GLvoid *vbo_data);
	void update_ebo(uint *ebo, int ebo_size, const GLvoid *ebo_data);

	// Rendering Vertex Buffer Object (VBO)
	void render_regular(bool dynamic = false);
	void render_wireframe(bool dynamic = false);
	void render_vertices(bool dynamic = false);
	void render_as_points(bool dynamic = false);
	void render_depth(bool dynamic = false);
	void render(bool dynamic = false);

	// State of VBO
	bool isDirty;
	void setDirty(bool state);
	bool isReady;
	bool isVBOEnabled;

	// Rendering flags
	bool isFlatShade;
	RENDER_MODE	render_mode;
	void setRenderMode(RENDER_MODE r);
};
