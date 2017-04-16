TEMPLATE = app
TARGET = Workspace
QT += core gui xml opengl
CONFIG += debug

INCLUDEPATH += ./GeneratedFiles \
./GeneratedFiles/Debug \
./GraphicsLibrary/Mesh/SurfaceMesh \
./Utility

include(Workspace.pri)

win32{
        LIBS += -L$$PWD"/../libQGLViewer/QGLViewer" \
        -lopengl32 \
        -lglu32

        CONFIG(debug, debug|release){
                LIBS += -lQGLViewerd2
        }else{
                LIBS += -lQGLViewer2
        }
}

linux {
    CONFIG(debug, debug|release){
            QMAKE_CXXFLAGS+= -ggdb -g3 -O0 -fopenmp
            QMAKE_LFLAGS *= -fopenmp
            LIBS += -lGLEW -lGLU -lGL -lQGLViewer
    }
}

macx{
    LIBS += -lGLEW
}

QT_VERSION=$$[QT_VERSION]

### Unix configuration ###
unix {
  CONFIG -= debug debug_and_release
  CONFIG *= release

  isEmpty( PREFIX ) {
    # Try same INCLUDE_DIR and LIB_DIR parameters than for the make install.
    PREFIX=/usr
  }

  # LIB_NAME
  LIB_NAME = libQGLViewer*.so*
  macx|darwin-g++ {
    LIB_NAME = libQGLViewer.$${QMAKE_EXTENSION_SHLIB}
  }
  hpux {
    LIB_NAME = libQGLViewer*.sl*
  }

  !isEmpty( QGLVIEWER_STATIC ) {
    LIB_NAME = libQGLViewer*.a
  }

  # LIB_DIR
  isEmpty( LIB_DIR ) {
    LIB_DIR = $${PREFIX}/lib
  }

  !exists( $${LIB_DIR}/$${LIB_NAME} ) {
    exists( ../../QGLViewer/$${LIB_NAME} ) {
      #message( The library was found in ../../QGLViewer which will be set as the LIB_DIR )
      LIB_DIR = ../../QGLViewer
    }
  }

  !exists( $${LIB_DIR}/$${LIB_NAME} ) {
    exists( ../../QGLViewer-build-desktop/$${LIB_NAME} ) {
      #message( The library was found in ../../QGLViewer-build-desktop which will be set as the LIB_DIR )
      LIB_DIR = ../../QGLViewer-build-desktop
    }
  }

  macx|darwin-g++ {
    !exists( $${LIB_DIR}/$${LIB_NAME} ) {
      # DYLIB was not found, try to find Framework instead
      LIB_NAME = QGLViewer.framework/QGLViewer
      LIB_DIR = ~/Library/Frameworks
      # qmake does not handle tilde
      LIB_DIR = $$system(cd $${LIB_DIR};pwd)

      !exists( $${LIB_DIR}/$${LIB_NAME} ) {
        exists( ../../QGLViewer/$${LIB_NAME} ) {
          #message( The framework was found in ../../QGLViewer which will be set as the LIB_DIR )
          LIB_DIR = ../../QGLViewer
        }
      }

      !exists( $${LIB_DIR}/$${LIB_NAME} ) {
        exists( ../../QGLViewer-build-desktop/$${LIB_NAME} ) {
          #message( The framework was found in ../../QGLViewer-build-desktop which will be set as the LIB_DIR )
          LIB_DIR = ../../QGLViewer-build-desktop
        }
      }
    }
  }

  !exists( $${LIB_DIR}/$${LIB_NAME} ) {
    message( Unable to find $${LIB_NAME} in $${LIB_DIR}. Make sure you have built it. )
    message( You should run qmake LIB_DIR=/path/to/QGLViewer/$${LIB_NAME} )
  }

  # The actual directory where the library/framework was found
  LIB_DIR_ABSOLUTE_PATH = $$system(cd $${LIB_DIR};pwd)

  # INCLUDE_DIR
  isEmpty( INCLUDE_DIR ) {
    INCLUDE_DIR = $${PREFIX}/include
  }

  macx|darwin-g++ {
          !exists( $${INCLUDE_DIR}/QGLViewer/qglviewer.h ) {
          INCLUDE_DIR=$${LIB_DIR}/QGLViewer.framework
              exists( $${LIB_DIR}/QGLViewer.framework/Headers/QGLViewer/qglviewer.h ) {
                 INCLUDE_DIR = $${LIB_DIR}/QGLViewer.framework/Headers
              }
          }
  }

  !exists( $${INCLUDE_DIR}/QGLViewer/qglviewer.h ) {
    exists( ../../QGLViewer/qglviewer.h ) {
      # message( libQGLViewer header files were not installed in standard $${INCLUDE_DIR} directory )
      # message( Headers were found in ../.. which will be set as the INCLUDE_DIR )
      INCLUDE_DIR = ../..
    }
  }

  !exists( $${INCLUDE_DIR}/QGLViewer/qglviewer.h ) {
    message( Unable to find QGLViewer/qglviewer.h in $${INCLUDE_DIR} )
    error( Use qmake INCLUDE_DIR=/path/to/QGLViewerDir )
  }


  macx|darwin-g++ {
    # On Mac, the lib path can be specified in the executable using install_name_tool
    contains( LIB_NAME, ".*QGLViewer.framework.*" ) {
      # If framework was not found in a standard directory
      !contains( LIB_DIR, ".*/Library/Frameworks/*$" ) {
        QMAKE_LFLAGS += -F$${LIB_DIR}
        !plugin:QMAKE_POST_LINK=install_name_tool -change QGLViewer.framework/Versions/2/QGLViewer $${LIB_DIR_ABSOLUTE_PATH}/QGLViewer.framework/Versions/2/QGLViewer $${TARGET}.app/Contents/MacOS/$${TARGET}
      }
      LIBS += -F$${LIB_DIR} -framework QGLViewer
    } else {
        !plugin:QMAKE_POST_LINK=install_name_tool -change libQGLViewer.2.dylib $${LIB_DIR_ABSOLUTE_PATH}/libQGLViewer.2.dylib $${TARGET}.app/Contents/MacOS/$${TARGET}
        LIBS *= -L$${LIB_DIR} -lQGLViewer
    }
  } else {
    isEmpty(QMAKE_LFLAGS_RPATH) {
      !plugin:QMAKE_LFLAGS += -Wl,-rpath,$${LIB_DIR_ABSOLUTE_PATH}
    } else {
      !plugin:QMAKE_RPATHDIR *= $${LIB_DIR_ABSOLUTE_PATH}
    }
    LIBS *= -L$${LIB_DIR} -lQGLViewer

        # Qt 4.8 removed the GLU dependency
    QMAKE_LIBS_OPENGL *= -lGLU
  }

  # Paths were correctly detected
  INCLUDEPATH *= $${INCLUDE_DIR}
  DEPENDPATH  *= $${INCLUDE_DIR}

  !isEmpty( QGLVIEWER_STATIC ) {
    LIBS *= $${LIB_DIR}/$${LIB_NAME}
  }

  macx|darwin-g++ {
    !contains( QT_VERSION, "^4.*" ) {
      # Qt3 only
      LIBS *= -lobjc
      CONFIG -= thread
    }
  }

  # Remove debugging options in release mode makes files much smaller
  release:QMAKE_CFLAGS_RELEASE -= -g
  release:QMAKE_CXXFLAGS_RELEASE -= -g

  # Intermediate files are created in an hidden folder
  MOC_DIR = .moc
  OBJECTS_DIR = .obj
}
