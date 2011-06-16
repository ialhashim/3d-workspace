#pragma once

#include "Mesh.h"

// Helper structs
struct SamplePoint{
	Vec pos, n;
	double weight;
	double u,v;
	int findex; // index of sampled face
	int flag;

	SamplePoint(const Vec& position = Vec(), const Vec& normal = Vec(), 
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
	Face * f;

	AreaFace(double a = 0.0, Face * face = NULL) : area(a), f(face){};

	bool operator< (const AreaFace & af) const { return area < af.area; }
	void setValue (double val) { area = val; }
};

enum SamplingMethod { FACE_CENTER, RANDOM_BARYCENTRIC };

// Class definition
class Sampler{

private:

public:
	
	SamplingMethod method;

	Sampler( Mesh * srcMesh = NULL, SamplingMethod samplingMethod = RANDOM_BARYCENTRIC );

	// Get samples
	SamplePoint getSample();
	Vector<SamplePoint> getSamples(int numberSamples);

	// Bias samples
	void resampleWithBias();
	void clearBias();
	Vector<double> bias;

	Mesh * mesh;
	double totalMeshArea;

	// For Monte Carlo
	Vector<AreaFace> interval;
	Vector<double> faceAreas;
	Vector<double> faceProbability;

	// DEBUG:
	void draw(const Vector<SamplePoint> & samples);

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
