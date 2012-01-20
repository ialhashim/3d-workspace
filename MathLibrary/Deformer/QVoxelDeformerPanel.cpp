#include "QVoxelDeformerPanel.h"

QVoxelDeformerPanel::QVoxelDeformerPanel()
{
	activeDeformer = NULL;
	activeScene = NULL;

	dw.setupUi(this);

	// Connections
	connect(dw.buildButton, SIGNAL(clicked()), SLOT(onBuildButtonClicked()));
	connect(dw.falloff, SIGNAL(valueChanged(double)), SLOT(onFalloffChanged(double)));
}

void QVoxelDeformerPanel::setActiveScene( Scene * newScene)
{
	activeScene = newScene;

	if(!activeScene) return;

	connect(this, SIGNAL(updateScene()), activeScene, SLOT(updateActiveObject()));
}

void QVoxelDeformerPanel::onBuildButtonClicked()
{
	if(!activeScene || !activeScene->activeObject())
		return;

	Primitive * prim = activeScene->activeObject()->controller->getSelectedPrimitive();

	if(!prim) return;

	QSurfaceMesh* mesh = prim->m_mesh;

	activeDeformer = new VoxelDeformer(mesh, dw.voxelSize->value());

	connect(activeDeformer, SIGNAL(meshDeformed()), activeScene, SLOT(updateActiveObject()));

	emit( deformerCreated(activeDeformer) );
}

void QVoxelDeformerPanel::onFalloffChanged(double value)
{
	if(!activeDeformer) return;

	activeDeformer->sigmaControl = value;

	emit(updateScene());
}
