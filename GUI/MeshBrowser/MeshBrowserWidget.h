#pragma once

#include <QStringList>
#include <QVector>
#include "ui_MeshBrowserForm.h"
#include "QuickMeshViewer.h"

class MeshBrowserWidget : public QDialog{
	Q_OBJECT

public:
	MeshBrowserWidget();

	QString selectedFile();

private:
	Ui::MeshBrowserDialog dialog;
	QWidget * thumbWidget;

	QVector< QVector<QuickMeshViewer*> > viewers;
	QuickMeshViewer* activeViewer;

	QString path;
	QStringList files;

	int countX, countY;
	int numActiveViewers;

protected:
	virtual void showEvent(QShowEvent * event);

public slots:
	void changePath();
	void loadMeshes(QString using_path);
	void loadMeshes();
	void showNumViewers(int n);
	void loadCurrentMeshes();
	void setActiveViewer(QuickMeshViewer*);
	void refresh();

signals:
	void pathChanged(QString);
};

#include <QThread>
class LoaderThread : public QThread{
	Q_OBJECT
public:
	LoaderThread(QuickMeshViewer *, QString);
	QuickMeshViewer * viewer;
	QString fileName;
protected:
	void run();
};
