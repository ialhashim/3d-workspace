#include "Sampler.h"

#include "SimpleDraw.h"
#include "Stats.h"

Sampler::Sampler(Mesh * srcMesh, SamplingMethod samplingMethod)
{
	isReady = false;

	if(srcMesh == NULL) 
		return;
	else
		mesh = srcMesh;

	method = samplingMethod;

	// Sample based on method selected
	if( method == RANDOM_BARYCENTRIC )
	{
		// Compute all faces area
		faceAreas = Vector<double> (mesh->numberOfFaces());
		faceProbability = Vector<double> (mesh->numberOfFaces());
		totalMeshArea = 0;

		for(StdList<Face>::iterator f = mesh->face.begin(); f != mesh->face.end(); f++)
		{
			faceAreas[f->index] = f->area();
			totalMeshArea += faceAreas[f->index];
		}

		for(int i = 0; i < (int) faceAreas.size(); i++)
			faceProbability[i] = faceAreas[i] / totalMeshArea;

		interval = Vector<AreaFace>(mesh->numberOfFaces() + 1);
		interval[0] = AreaFace(0.0, mesh->faceIndexMap[0]);
		int i = 0;

		// Compute mesh area in a cumulative manner
		for(int j = 0; j < (int) faceAreas.size(); j++)
		{
			interval[i+1] = AreaFace(interval[i].area + faceProbability[j], mesh->faceIndexMap[j]);
			i++;
		}

		// For importance sampling
		clearBias();
	}
	else if( method ==  FACE_CENTER )
	{
		// No preparations needed..
	}

	isReady = true;
}

SamplePoint Sampler::getSample()
{
	SamplePoint sp;

	if( method == RANDOM_BARYCENTRIC )
	{
		// r, random point in the area
		double r = uniform();

		// Find corresponding face
		Vector<AreaFace>::iterator it = lower_bound(interval.begin(), interval.end(), AreaFace(Min(r,interval.back().area)));
		Face * face = it->f;

		// Add sample from that face
		double b[3]; RandomBaricentric(b);

		sp = SamplePoint( face->getBary(b[0], b[1]), *mesh->fn(face->index), 1.0, face->index, b[0], b[1]);
	}
	else if( method ==  FACE_CENTER )
	{
		int fcount = mesh->numberOfFaces();

		int randTriIndex = (int) (fcount * (((double)rand()) / (double)RAND_MAX)) ;

		if( randTriIndex >= fcount )
			randTriIndex = fcount - 1;

		Face * tri = mesh->f( randTriIndex );

		// Get triangle center and normal
		Vec triCenter = tri->center();
		Vec triNor = *mesh->fn(randTriIndex);

		sp = SamplePoint(triCenter, triNor, tri->area(), tri->index);
	}

	return sp;
}

Vector<SamplePoint> Sampler::getSamples(int numberSamples)
{
	Vector<SamplePoint> samples;

	for(int i = 0; i < numberSamples; i++)
	{
		samples.push_back( getSample() );
	}

	return samples;
}

void Sampler::clearBias()
{
	bias = faceProbability;
}

void Sampler::resampleWithBias()
{
	isReady = false;

	// Fill with lowest bias
	double minBias = *min_element(bias.begin(), bias.end());
	for(int i = 0; i < (int) bias.size(); i++)	
		bias[i] = Max(bias[i], minBias);

	double totalNewArea = 0;

	for(int i = 0; i < (int) faceAreas.size(); i++)
	{
		faceAreas[i] = bias[i];
		totalNewArea += faceAreas[i];
	}

	for(int i = 0; i < (int) faceAreas.size(); i++)
		faceProbability[i] = faceAreas[i] / totalNewArea;

	interval = Vector<AreaFace>(mesh->numberOfFaces() + 1);
	interval[0] = AreaFace(0.0, mesh->faceIndexMap[0]);
	int i = 0;

	// Compute mesh area in a cumulative manner
	for(int j = 0; j < (int) faceAreas.size(); j++)
	{
		interval[i+1] = AreaFace(interval[i].area + faceProbability[j], mesh->faceIndexMap[j]);
		i++;
	}

	// Normalize new total mesh area
	double total = interval.back().area;
	for(int i = 0; i < (int)interval.size(); i++)
		interval[i].area /= total;

	isReady = true;
}

void Sampler::draw(const Vector<SamplePoint> & samples)
{
	glEnable(GL_BLEND); 
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glPointSize(4.0);

	glBegin(GL_POINTS);
	for(int i = 0; i < (int)samples.size(); i++)
	{
		const SamplePoint * p = &samples[i];

		glColor4f(1,0,0, p->weight);
		glNormal3fv(p->n);
		glVertex3dv(p->pos);
	}
	glEnd();
}
