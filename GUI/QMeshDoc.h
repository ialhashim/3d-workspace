#pragma once

#include <QObject>
#include <QMap>
#include "GraphicsLibrary/Mesh/QSegMesh.h"

class Workspace;

class QMeshDoc : public QObject
{
	Q_OBJECT

public:
	QMeshDoc();
	~QMeshDoc();

	// Objects
	QMap<QString, QSegMesh*> all_objects;
	uint global_id;
	QSegMesh * getObject( QString objectId );

	QString loadFile(QString fileName);

	Workspace * workspace;

signals:
	void objectImported( QSegMesh* mesh );
	void printMessage( QString message );

public slots:
	void importObject( );
	void exportObject( QSegMesh * mesh );
	void deleteObject( QString objectId );
};

