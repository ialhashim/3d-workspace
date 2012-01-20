#pragma once

#include "QSurfaceMesh.h"
#include "QSegMesh.h"

// Helper structures
struct SamplePoint{
	Vec3d pos, n;
	double weight;
	double u,v;
	int findex; // index of sampled face
	int flag;

	SamplePoint(const Vec3d& position = Vec3d(), const Vec3d& normal = Vec3d(), 
		double Weight = 0.0, int face_index = -1.0, double U = 0.0, double V = 0.0, int flags = 0)
	{
		pos = position;
		n = normal;
		weight = Weight;
		findex = face_index;
		u = U;
		v = V;
		flag = flags;
	}
};

struct AreaFace{
	double area;
	Surface_mesh::Face f;

	AreaFace(double a = 0.0, Surface_mesh::Face face = Surface_mesh::Face()) : area(a), f(face){};

	bool operator< (const AreaFace & af) const { return area < af.area; }
	void setValue (double val) { area = val; }
};

enum SamplingMethod { FACE_CENTER, RANDOM_BARYCENTRIC };

// Class definition
class Sampler{

private:

public:
	
	SamplingMethod method;

	Sampler(QSurfaceMesh * srcMesh = NULL, SamplingMethod samplingMethod = RANDOM_BARYCENTRIC );
	Sampler(void * srcMesh, SamplingMethod samplingMethod);
	// Get samples
	SamplePoint getSample(double weight = 0.0);
	StdVector<SamplePoint> getSamples(int numberSamples, double weight = 0.0);

	static StdVector<SamplePoint> getSamplesFromQSegMesh(QSegMesh* srcMesh, int numberSamples);
	// Bias samples
	void resampleWithBias();
	void clearBias();
	StdVector<double> bias;

	QSurfaceMesh * mesh;
	double totalMeshArea;

	// For Monte Carlo
	StdVector<AreaFace> interval;
	StdVector<double> faceAreas;
	StdVector<double> faceProbability;

	// DEBUG:
	static void draw(const StdVector<SamplePoint> & samples);

	bool isReady;
};

// Helper functions
static inline void RandomBaricentric(double * interp)
{
	interp[1] = uniform();
	interp[2] = uniform();

	if(interp[1] + interp[2] > 1.0)
	{
		interp[1] = 1.0 - interp[1];
		interp[2] = 1.0 - interp[2];
	}

	interp[0] = 1.0 - (interp[1] + interp[2]);
}
