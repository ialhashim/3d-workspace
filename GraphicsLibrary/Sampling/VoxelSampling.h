#pragma once

#include "GraphicsLibrary/Mesh/QSegMesh.h"
#include "GraphicsLibrary/Basic/Plane.h"
#include "GraphicsLibrary/Voxel/Voxeler.h"

class VoxelSampling{
public:
	static std::vector<Vec3d> sample(QSurfaceMesh * m, double voxel_size)
	{
		return Voxeler (m, voxel_size).getVoxelCenters();
	}
};
