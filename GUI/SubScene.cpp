#include "SubScene.h"
#include "SimpleDraw.h"

SubScene::SubScene( Scene * parentScene, int X, int Y, int Width, int Height )
{
	activeScene = parentScene;

	x = X, y = Y;
	width = Width; height = Height;

	if(parentScene != NULL)	isDraw = true;
	else isDraw = false;

	isSelected = false;

	caption = "SubScene";
}

QSegMesh * SubScene::activeObject()
{
	if(activeScene)
		return activeScene->activeObject();
	else
		return NULL;
}

void SubScene::drawFrame(Vec4d color)
{
	activeScene->startScreenCoordinatesSystem();
	
	double z = 0.00001;

	Vec3d v1(x,y,z);
	Vec3d v2(x+ width,y,z);
	Vec3d v3(x+ width,y+height,z);
	Vec3d v4(x,y+height,z);

	SimpleDraw::DrawSquare(v1,v2,v3,v4, true,3,color[0],color[1],color[2], color[3]);
	
	activeScene->stopScreenCoordinatesSystem();
}

void SubScene::draw()
{
	if(!isDraw) return;

	// Color
	Vec4d selectedColor(1,1,0,0.2), regularColor(0.8,0.8,0.9, 0.1);
	Vec4d c = regularColor;
	
	if(isSelected)
		c = selectedColor;

	drawFrame(c);
}

bool SubScene::contains( Vec2i pixel )
{
	return RANGE(pixel.x(), x, x + width) && RANGE(pixel.y(), y, y + height);
}

void SubScene::postDraw()
{
	glColor4dv(Vec4d(0,0,0,0.4));
	activeScene->renderText(x + 6, y + 16, caption);

	glColor4dv(Vec4d(1,1,1,1));
	activeScene->renderText(x + 5, y + 15, caption);
}
