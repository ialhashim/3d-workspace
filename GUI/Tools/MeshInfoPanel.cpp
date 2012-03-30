#include "MeshInfoPanel.h"

MeshInfoPanel::MeshInfoPanel()
{
	infoWidget.setupUi(this);
}

void MeshInfoPanel::setActiveScene( Scene * scene )
{
	if(activeScene != scene)
		activeScene = scene;

	if(!activeScene || !activeObject())
	{
		QString zero = QString::number(0);

		infoWidget.numVertices->setText(zero);
		infoWidget.numFaces->setText(zero);
		infoWidget.numSegments->setText(zero);

		return;
	}
		
	infoWidget.numVertices->setText(QString::number(activeObject()->nbVertices()));
	infoWidget.numFaces->setText(QString::number(activeObject()->nbFaces()));
	infoWidget.numSegments->setText(QString::number(activeObject()->nbSegments()));
}

QSegMesh* MeshInfoPanel::activeObject()
{
	if (activeScene)
		return activeScene->activeObject();
	else 
		return NULL;
}
