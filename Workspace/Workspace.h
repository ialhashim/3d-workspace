#ifndef WORKSPACE_H
#define WORKSPACE_H

#include <QtGui/QMainWindow>
#include "ui_Workspace.h"

#include <QMdiSubWindow>
#include <QFileDialog>

class Workspace : public QMainWindow
{
	Q_OBJECT

public:
	Workspace(QWidget *parent = 0, Qt::WFlags flags = 0);
	~Workspace();

public slots:
	void addNewScene();
	void importObject();

signals:
	void importedObject(QString fileName);

private:
	Ui::WorkspaceClass ui;
};

#endif // WORKSPACE_H
