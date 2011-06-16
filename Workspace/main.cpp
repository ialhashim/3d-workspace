#include "Workspace.h"
#include <QtGui/QApplication>
#include <QGLFormat>

double Epsilon = 1.0e-7f;

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);

	// Anti-aliasing
	QGLFormat glf = QGLFormat::defaultFormat();
	glf.setSamples(8);
	QGLFormat::setDefaultFormat(glf);

	Workspace w;
	w.show();
	return a.exec();
}
