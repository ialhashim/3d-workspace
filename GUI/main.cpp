#include "global.h"
#include "Workspace.h"
#include <QApplication>
#include <QGLFormat>
#include <QApplication>
#include <QDesktopWidget>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);

	DEFAULT_FILE_PATH = "";

	// Anti-aliasing
	QGLFormat glf = QGLFormat::defaultFormat();
	glf.setSamples(8);
	QGLFormat::setDefaultFormat(glf);

	// Create main window
	Workspace w;
	w.move(QApplication::desktop()->availableGeometry().center() - w.rect().center());
	w.show();
	w.addNewScene();

	return a.exec();
}
