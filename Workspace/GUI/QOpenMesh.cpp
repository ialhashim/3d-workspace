#include "QOpenMesh.h"

// OpenGL stuff
#include "OpenMeshGL.h"

void QOpenMesh::init()
{
	color = QColor(255,255,255,255);

	computeNormals();

	computeBoundingBox();
}

void QOpenMesh::computeNormals()
{
	// Add properties to our mesh container
	request_face_normals();
	request_vertex_normals();

	// Compute the normals per face, then per vertex
	update_normals();
}

void QOpenMesh::computeBoundingBox()
{
	ConstVertexIter v_it(vertices_begin()), v_end (vertices_end());

	Point firstPoint = point(v_it);

	bbox_min = bbox_max = firstPoint;

	for (; v_it != v_end; ++v_it )
	{               
		Point p = point(v_it);

		for (int i = 0; i < 3; i++){
			bbox_min[i] = ((bbox_min[i]<p[i])?bbox_min[i]:p[i]);
			bbox_max[i] = ((bbox_max[i]>p[i])?bbox_max[i]:p[i]);
		}
	}
}

void QOpenMesh::draw()
{
	ConstFaceIter f_it(faces_sbegin()), f_end(faces_end());
	ConstFaceVertexIter v_it;

	glEnable(GL_LIGHTING);

	// Entire mesh color (for now)
	glColor4f(color.redF(), color.greenF(), color.blueF(), color.alphaF());

	glBegin(GL_TRIANGLES);
	for (; f_it != f_end; ++f_it)
	{
		v_it = cfv_iter(f_it);

		glNormal(normal(++v_it));
		glVertex(point(v_it));

		glNormal(normal(++v_it));
		glVertex(point(v_it));

		glNormal(normal(++v_it));
		glVertex(point(v_it));
	}
	glEnd();

}

double QOpenMesh::radius()
{
	return 0.5 * (bbox_max - bbox_min).length();
}
