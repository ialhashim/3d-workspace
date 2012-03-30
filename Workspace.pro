TEMPLATE = app
TARGET = Workspace
DESTDIR = ./
QT += core gui xml opengl
CONFIG += debug console

DEFINES += QT_XML_LIB QT_OPENGL_LIB qh_QHpointer QT_DLL

INCLUDEPATH += ./GeneratedFiles \
    ./GeneratedFiles/Debug \
    ./GraphicsLibrary/Mesh/SurfaceMesh \
    ./Utility \


win32{
	LIBS += -L$$PWD"/GUI/Viewer/libQGLViewer/QGLViewer/lib" \
	    -lopengl32 \
	    -lglu32 \
	    -lQGLViewerd2
}

unix {
	CONFIG(debug, debug|release){
	    QMAKE_CXXFLAGS+= -ggdb -g3 -O0 -fopenmp
	    QMAKE_LFLAGS *= -fopenmp
	    LIBS += -lGLEW -lGLU -lGL -lQGLViewer
	}
}


INCLUDEPATH += ./GeneratedFiles \
    ./GeneratedFiles/Debug \
    . \
    ./GraphicsLibrary/Mesh/SurfaceMesh \
    ./Utility

DEPENDPATH += .
MOC_DIR += ./GeneratedFiles/debug
OBJECTS_DIR += debug
UI_DIR += ./GeneratedFiles
RCC_DIR += ./GeneratedFiles
include(Workspace.pri)
