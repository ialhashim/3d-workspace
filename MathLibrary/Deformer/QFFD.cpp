#include "QFFD.h"
#include "Utility/SimpleDraw.h"

FFD current_ffd;

QFFD::QFFD( QSurfaceMesh * src_mesh, FFD_FitType fit_type, Vec3i gridResolution )
{
	isReady = false;

	m_mesh = src_mesh;

	current_ffd = FFD(src_mesh, fit_type);

	// Fit control points using bounding box
	ffd()->bbFit( gridResolution );

	// Connect control points to update function
	foreach(QControlPoint * cp, this->ffd()->points)
	{
		connect(cp, SIGNAL(manipulated()), this, SLOT(updateMesh()));
	}

	isReady = true;
}

void QFFD::draw()
{
	if(!isReady) return;

	// Visual representation
	foreach(QControlPoint * cp, this->ffd()->points)
	{
		SimpleDraw::IdentifyPoint(cp->pos, 1.0f, 0.6f, 0, 8.0f);
		
		drawConnections(cp);
	}

	// Draw sphere around selected control points
	foreach(int idx, selectedPoints){
		glColor4dv(Color(0.7,0.9,0,1)); 
		SimpleDraw::DrawSphere(ffd()->points[idx]->pos, Min(ffd()->width, Min(ffd()->length, ffd()->width)) * 0.02);
	}

	// Debug
	foreach(Vec3d db_point, this->ffd()->dbPoints)
		SimpleDraw::IdentifyPoint(db_point, 0.6f, 0.5f, 0.1f, 4.0f);

	for(StdVector< Pair<Vec3d,Vec3d> >::iterator db_line = ffd()->dbLines.begin(); 
		db_line != ffd()->dbLines.end(); db_line++ )
	{
		SimpleDraw::IdentifyLine(db_line->first, db_line->second);
	}
}

void QFFD::drawNames()
{
	glPointSize(20);

	foreach(QControlPoint * cp, this->ffd()->points)
	{
		glPushName(cp->idx);
		glBegin(GL_POINTS);
		glVertex3dv(cp->pos);
		glEnd();
		glPopName();
	}
}

void QFFD::drawConnections(QControlPoint * cp)
{
	Vec3d p = cp->pos, q;

	// Index on grid
	Vec3i idx = cp->gridIdx;
	int x = idx.x(), y = idx.y(), z = idx.z();

	// Set line width
	glLineWidth(2.0);

	// Set color
	glEnable(GL_BLEND); 
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glColor4dv(Color(0.992, 0.584, 0, 0.6)); // orange-ish

	glDisable(GL_LIGHTING);

	glBegin(GL_LINES);
	if(x > 0){
		q = ffd()->points[ffd()->pointsGridIdx[x - 1][y][z]]->pos;
		glVertex3dv(p);	glVertex3dv(q);
	}

	if(y > 0){
		q = ffd()->points[ffd()->pointsGridIdx[x][y - 1][z]]->pos;
		glVertex3dv(p);	glVertex3dv(q);
	}

	if(z > 0){
		q = ffd()->points[ffd()->pointsGridIdx[x][y][z - 1]]->pos;
		glVertex3dv(p);	glVertex3dv(q);
	}
	glEnd();

	glEnable(GL_LIGHTING);
}

void QFFD::postSelection(int idx)
{
	if(idx == -1){
		selectedPoints.clear();
		return;
	}

	selectedPoints.push_back(idx);
}

QControlPoint * QFFD::getQControlPoint( int index )
{
	return ffd()->points[index];
}

FFD * QFFD::ffd()
{
	return &current_ffd;
}

void QFFD::updateMesh()
{
	ffd()->apply();

	emit(meshDeformed());
}
