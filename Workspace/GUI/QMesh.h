#pragma once
#include <QObject>

#undef __gl_h_
#include "Mesh.h"

class QMesh : public QObject, public Mesh
{
	Q_OBJECT

public:
	QMesh() : Mesh(){ }
};
