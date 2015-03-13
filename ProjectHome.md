A simple GUI testbed for 3D processing and visualization. Cross-platform (Windows, OSX, Linux). Source code is written in C++. Required library for 3D viewer http://www.libqglviewer.com/.

(Note that there are **no GUI** tools ready for you to use these. You have to write your own as in the rotation example provided.)

Implemented algorithms:
  * **Mesh simplification** using quadric error
  * **Sampling** (regular recursive, Monte Carlo, sphere packing, voxels)
  * **Skeleton extraction**
  * **Smoothing** (Laplacian, scale-dependent, mean curvature flow)
  * **Space partitioning** (octree, kdtree)
  * **Mesh voxelization**
  * **Subdivision** (loop, modified butterfly, longest edge)
  * **Minimum oriented bounding box (OBB)**
  * **Coordinates** (mean value coordinates, Green coordinates)
  * **Curvature** (polynomial fitting, two other implementations)
  * **Deformation** (As-rigid-as-possible ([ARAP](http://www.youtube.com/watch?v=owTuTyz41mM)), [Linear Rotation-Invariant Coordinates](http://www.youtube.com/watch?v=bdRzRW5G3HA), FFD, voxel deformation)
  * **Skinning** with dual quaternions
  * **Colormap**
  * **Mesh Browser** display 3D thumbnails of all meshes in a folder
  * [re-meshing](http://dl.acm.org/citation.cfm?id=1057457)

![http://www.cs.sfu.ca/~iaa7/personal/3d-workspace.png](http://www.cs.sfu.ca/~iaa7/personal/3d-workspace.png)