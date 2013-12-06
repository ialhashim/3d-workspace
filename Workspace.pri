HEADERS += ./MathLibrary/Bounding/BoundingBox.h \
    ./Utility/ColorMap.h \
    ./GraphicsLibrary/Sampling/EdgeSampler.h \
    ./Utility/Graph.h \
    ./Utility/HashTable.h \
    ./GraphicsLibrary/Mesh/SurfaceMesh/IO_.h \
    ./GraphicsLibrary/Basic/Line.h \
    ./Utility/Macros.h \
    ./Utility/NoiseGen.h \
    ./GraphicsLibrary/SpacePartition/Octree.h \
    ./GraphicsLibrary/Basic/Plane.h \
    ./GraphicsLibrary/Basic/PolygonArea.h \
    ./GraphicsLibrary/Mesh/SurfaceMesh/Quadric.h \
    ./GraphicsLibrary/Sampling/RegularRecursive.h \
    ./GraphicsLibrary/Sampling/Sampler.h \
    ./Utility/SimpleDraw.h \
    ./Utility/Sleeper.h \
    ./GraphicsLibrary/Sampling/SliceSampling.h \
    ./GraphicsLibrary/Sampling/SpherePackSampling.h \
    ./Utility/Stats.h \
    ./GraphicsLibrary/Mesh/SurfaceMesh/Surface_mesh.h \
    ./GraphicsLibrary/Basic/Triangle.h \
    ./GL/VBO/VBO.h \
    ./GraphicsLibrary/Mesh/SurfaceMesh/Vector.h \
    ./GraphicsLibrary/Voxel/Voxel.h \
    ./GraphicsLibrary/Sampling/VoxelSampling.h \
    ./GraphicsLibrary/Voxel/Voxeler.h \
    ./GUI/global.h \
    ./GraphicsLibrary/Mesh/SurfaceMesh/properties.h \
    ./GUI/Tools/MeshInfoPanel.h \
    ./GUI/QMeshDoc.h \
    ./GraphicsLibrary/Mesh/QSegMesh.h \
    ./GraphicsLibrary/Mesh/QSurfaceMesh.h \
    ./GUI/Scene.h \
    ./GUI/Tools/TransformationPanel.h \
    ./GUI/Workspace.h \
    GUI/MeshBrowser/QuickMeshViewer.h \
    GUI/MeshBrowser/QuickMesh.h \
    GUI/MeshBrowser/MeshBrowserWidget.h \
    GL/GLee.h
SOURCES += ./MathLibrary/Bounding/BoundingBox.cpp \
    ./Utility/ColorMap.cpp \
    ./GraphicsLibrary/Mesh/SurfaceMesh/IO_.cpp \
    ./GraphicsLibrary/Mesh/SurfaceMesh/IO_obj.cpp \
    ./GraphicsLibrary/Mesh/SurfaceMesh/IO_off.cpp \
    ./GraphicsLibrary/Mesh/SurfaceMesh/IO_stl.cpp \
    ./GraphicsLibrary/Basic/Line.cpp \
    ./GUI/Tools/MeshInfoPanel.cpp \
    ./GraphicsLibrary/SpacePartition/Octree.cpp \
    ./GraphicsLibrary/Basic/Plane.cpp \
    ./GUI/QMeshDoc.cpp \
    ./GraphicsLibrary/Mesh/QSegMesh.cpp \
    ./GraphicsLibrary/Mesh/QSurfaceMesh.cpp \
    ./GraphicsLibrary/Sampling/Sampler.cpp \
    ./GUI/Scene.cpp \
    ./Utility/SimpleDraw.cpp \
    ./Utility/Stats.cpp \
    ./GraphicsLibrary/Mesh/SurfaceMesh/Surface_mesh.cpp \
    ./GUI/Tools/TransformationPanel.cpp \
    ./GraphicsLibrary/Basic/Triangle.cpp \
    ./GL/VBO/VBO.cpp \
    ./GraphicsLibrary/Voxel/Voxeler.cpp \
    ./GUI/Workspace.cpp \
    ./GUI/global.cpp \
    ./GraphicsLibrary/SpacePartition/kdtree.cpp \
    ./GUI/main.cpp \
    GUI/MeshBrowser/QuickMeshViewer.cpp \
    GUI/MeshBrowser/MeshBrowserWidget.cpp \
    GL/GLee.c
FORMS += ./GUI/Tools/MeshInfo.ui \
    ./GUI/Tools/RotationWidget.ui \
    ./GUI/Workspace.ui \
    GUI/MeshBrowser/MeshBrowserForm.ui
RESOURCES += Workspace.qrc
