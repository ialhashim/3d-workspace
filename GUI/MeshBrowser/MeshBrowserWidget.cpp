#include "Utility/Macros.h"
#include "MeshBrowserWidget.h"
#include "QuickMeshViewer.h"
#include <QFileDialog>

MeshBrowserWidget::MeshBrowserWidget()
{
	// Setup dialog from designed form
	dialog.setupUi(this);
	thumbWidget = dialog.thumbnailWidget;

	countX = 3;
	countY = 4;

	viewers.resize(countX);
	for(int i = 0; i < countX; i++){
		viewers[i].resize(countY);
		for(int j = 0; j < countY; j++){
			viewers[i][j] = new QuickMeshViewer;
			connect(viewers[i][j], SIGNAL(gotFocus(QuickMeshViewer*)), SLOT(setActiveViewer(QuickMeshViewer*)));
		}
	}

	numActiveViewers = 0;
	activeViewer = viewers[0][0];

	for(int i = 0; i < countX; i++)
		for(int j = 0; j < countY; j++)
			dialog.thumbLayout->addWidget(viewers[i][j], i, j);

	// Default path
	path = "";

	// Connections
	connect(dialog.pathButton, SIGNAL(clicked()), SLOT(changePath()));
	connect(this, SIGNAL(pathChanged(QString)), dialog.folderLabel, SLOT(setText(QString)));
	connect(this, SIGNAL(pathChanged(QString)), SLOT(loadMeshes(QString)));
	connect(dialog.scrollBar, SIGNAL(valueChanged(int)), SLOT(loadMeshes()));

	changePath();
}

void MeshBrowserWidget::changePath()
{
	path = QFileDialog::getExistingDirectory(this, "Select folder", path, QFileDialog::ShowDirsOnly);

	emit(pathChanged(path));
}

void MeshBrowserWidget::showEvent( QShowEvent * event )
{
	activeViewer->setFocus();
}

void MeshBrowserWidget::showNumViewers( int n )
{
	int count = 0, activeCount = 0;

	for(int i = 0; i < countX; i++){
		for(int j = 0; j < countY; j++){
			if(count++ < n)
			{
				viewers[i][j]->isActive = true;
				viewers[i][j]->clearMesh();
				viewers[i][j]->resetView();
				
				activeCount++;
			}
			else
				viewers[i][j]->isActive = false;
		}
	}

	refresh();

	numActiveViewers = activeCount;
}

void MeshBrowserWidget::loadMeshes(QString using_path)
{
	path = using_path;

	// Get list of files
	QStringList filters;
	filters << "*.obj" << "*.off";
	files = QDir(path).entryList(filters);

	int numPages = ceil(double(files.size()) / double(countX * countY)) - 1;

	dialog.scrollBar->setRange(0, numPages);
	dialog.scrollBar->setValue(0);

	if(files.size())
		loadMeshes();
	else
	{
		showNumViewers(0);
		refresh();
	}
}

void MeshBrowserWidget::loadMeshes()
{
	loadCurrentMeshes();
}

void MeshBrowserWidget::refresh()
{
	for(int i = 0; i < countX; i++)
		for(int j = 0; j < countY; j++)
			viewers[i][j]->updateGL();
}

void MeshBrowserWidget::loadCurrentMeshes()
{
	int numViewers = countX * countY;
	int index = dialog.scrollBar->value() * numViewers;
	int curActive = Min(numViewers, files.size() - index);

	showNumViewers(curActive);

	int c = 0;

	for(int i = 0; i < countX; i++){
		for(int j = 0; j < countY; j++)
		{
			viewers[i][j]->clearMesh();
			viewers[i][j]->resetView();
		}
	}

	for(int i = 0; i < countX; i++){
		for(int j = 0; j < countY; j++)
		{
			if(index + c > files.size() - 1) return;

			QString fileName = path + "\\" + files[index + c];

			new LoaderThread(viewers[i][j], fileName);

			c++; if(c > curActive) return;
		}
	}
}

void MeshBrowserWidget::setActiveViewer( QuickMeshViewer* v)
{
	activeViewer = v;
}

QString MeshBrowserWidget::selectedFile()
{
	if(activeViewer) return activeViewer->meshFileName();
	return "";
}

LoaderThread::LoaderThread(QuickMeshViewer * v, QString file_name)
{
	this->viewer = v;
	this->fileName = file_name;
	this->start();
}

void LoaderThread::run()
{
	this->viewer->loadMesh(fileName);
	QThread::exit();
}
