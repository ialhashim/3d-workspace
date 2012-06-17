#pragma once

#include <QString>
#include <QVector3D>
#include <QFile>
#include <QTextStream>
#include <qgl.h>
#include <float.h>

class QuickMesh : public QObject{
	Q_OBJECT

public:
	QuickMesh()
	{
		verts.clear();
		tris.clear();
		isLoading = true;
		fileName = "";
	}

	bool isLoading;

	QString fileName;

	void draw()
	{
		if(isLoading) return;

		glEnable(GL_LIGHTING);
		glColor3d(1,1,1);

		if(tris.size())
		{
			glBegin(GL_TRIANGLES);
			foreach(const QVector<int> tri, tris){
				if(tri[0] < 0 || tri[1] < 0 || tri[2] < 0) continue;

				QVector3D v1 = verts[tri[0]], v2 =verts[tri[1]], v3 = verts[tri[2]];
				QVector3D n = QVector3D::crossProduct((v2-v1).normalized(), (v3-v1).normalized()).normalized();

				glNormal3d(n.x(), n.y(), n.z());

				glVertex3d(v1.x(), v1.y(), v1.z());
				glVertex3d(v2.x(), v2.y(), v2.z());
				glVertex3d(v3.x(), v3.y(), v3.z());
			}
			glEnd();
		}
		else{
			// Point cloud
			glPointSize(2.0f);
			glDisable(GL_LIGHTING);
			glBegin(GL_POINTS);
			foreach(const QVector3D v, verts) glVertex3d(v.x(), v.y(), v.z());
			glEnd();
		}

		glDisable(GL_LIGHTING);
	}

public slots:

	void load(QString fileame)
	{
		fileName = fileame;

		clear();

		QString ext = fileName.right(3).toLower();

		if(ext == "obj") loadOBJ(fileName);
		if(ext == "off") loadOFF(fileName);

		postProcess();

		isLoading = false;
	}

	void clear()
	{
		verts.clear();
		tris.clear();

		isLoading = true;
	}

private:

	void loadOBJ(QString fileName)
	{
		QFile file(fileName); file.open(QIODevice::ReadOnly | QIODevice::Text);
		QTextStream in(&file);

		while (!in.atEnd()){
			QStringList v = in.readLine().split(" ", QString::SkipEmptyParts);

			if(v.size() < 4) continue;

			if(v[0] == "v")
			{
				verts.push_back(QVector3D(v[1].toDouble(),v[2].toDouble(),v[3].toDouble()));
			} 
			
			if(v[0] == "f")
			{
				tris.push_back(QVector<int>(3));
				tris.back()[0] = ( v[1].replace("/", " ").split(" ", QString::SkipEmptyParts)[0].toInt() - 1 );
				tris.back()[1] = ( v[2].replace("/", " ").split(" ", QString::SkipEmptyParts)[0].toInt() - 1 );
				tris.back()[2] = ( v[3].replace("/", " ").split(" ", QString::SkipEmptyParts)[0].toInt() - 1 );
			}
		}
	}

	void loadOFF(QString fileName)
	{
		QFile file(fileName); file.open(QIODevice::ReadOnly | QIODevice::Text);
		QTextStream in(&file);

		// skip first line
		in.readLine(); 
		
		// Read number of verts and tris
		QStringList lineOne = in.readLine().split(" ", QString::SkipEmptyParts);
		if(lineOne.size() != 3) return;
		int num_v = lineOne[0].toInt();
		int num_f = lineOne[1].toInt();

		// Read vertices
		for(int vi = 0; vi < num_v; vi++){
			QStringList v = in.readLine().split(" ", QString::SkipEmptyParts);
			if(v.size() == 3) verts.push_back(QVector3D(v[0].toDouble(),v[1].toDouble(),v[2].toDouble()));
			else break;
		}

		if(verts.size() != num_v) return;

		// Read faces
		for(int fi = 0; fi < num_f; fi++){
			if(in.atEnd()) return;
			QStringList f = in.readLine().split(" ", QString::SkipEmptyParts);

			if(f.size() < 4) continue;
			tris.push_back(QVector<int>());

			for(int i = 1; i < f.size(); i++) tris.back().push_back(f[i].toInt());
		}
	}

	void postProcess()
	{
		// compute bounding box
		QVector3D bbmin (DBL_MAX-1, DBL_MAX-1, DBL_MAX-1);
		QVector3D bbmax = -bbmin;
		QVector3D center(0,0,0);
		foreach(const QVector3D v, verts)
		{
			if(v.x() < bbmin.x()) bbmin.setX(v.x());
			if(v.y() < bbmin.y()) bbmin.setY(v.y());
			if(v.z() < bbmin.z()) bbmin.setZ(v.z());

			if(v.x() > bbmax.x()) bbmax.setX(v.x());
			if(v.y() > bbmax.y()) bbmax.setY(v.y());
			if(v.z() > bbmax.z()) bbmax.setZ(v.z());
		}
		center = (bbmin + bbmax) * 0.5;

		// Normalize and move to center
		QVector3D d = bbmax - bbmin;
		double s = (d.x() > d.y())? d.x():d.y();
		s = (s>d.z())? s: d.z();
		for(int vi = 0; vi < verts.size(); vi++) 
			verts[vi] = (verts[vi] - center) / s;
	}
	
	QVector< QVector3D > verts;
	QVector< QVector<int> > tris;
};
