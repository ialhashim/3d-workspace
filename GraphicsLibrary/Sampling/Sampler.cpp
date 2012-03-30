#include "Sampler.h"

#include "Utility/SimpleDraw.h"
#include "Utility/Stats.h"
#include "Utility/ColorMap.h"


Sampler::Sampler(QSurfaceMesh * srcMesh, SamplingMethod samplingMethod)
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
		faceAreas = StdVector<double> (mesh->n_faces());
		faceProbability = StdVector<double> (mesh->n_faces());
		this->totalMeshArea = 0;

		Surface_mesh::Face_iterator fit, fend = mesh->faces_end();
		for (fit = mesh->faces_begin(); fit != fend; ++fit)
		{
			uint f_id = Surface_mesh::Face(fit).idx();

			faceAreas[f_id] = mesh->faceArea(fit);
			totalMeshArea += faceAreas[f_id];
		}

		for(int i = 0; i < (int) faceAreas.size(); i++)
			faceProbability[i] = faceAreas[i] / totalMeshArea;

		interval = StdVector<AreaFace>(mesh->n_faces() + 1);
		interval[0] = AreaFace(0.0, mesh->getFace(0));
		int i = 0;

		// Compute mesh area in a cumulative manner
		for(int j = 0; j < (int) faceAreas.size(); j++)
		{
			interval[i+1] = AreaFace(interval[i].area + faceProbability[j], mesh->getFace(j));
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

SamplePoint Sampler::getSample(double weight)
{
	SamplePoint sp;
	double r;
	double b[3];

	if( method == RANDOM_BARYCENTRIC )
	{
		// r, random point in the area
		r = uniform();

		// Find corresponding face
		StdVector<AreaFace>::iterator it = lower_bound(interval.begin(), interval.end(), AreaFace(Min(r,interval.back().area)));
		Surface_mesh::Face face = it->f;
		uint face_id = face.idx();

		// Add sample from that face
		RandomBaricentric(b);

		sp = SamplePoint( mesh->getBaryFace(face, b[0], b[1]), mesh->fn(face), weight, face_id, b[0], b[1]);
	}
	else if( method ==  FACE_CENTER )
	{
		int fcount = mesh->n_faces();

		int randTriIndex = (int) (fcount * (((double)rand()) / (double)RAND_MAX)) ;

		if( randTriIndex >= fcount )
			randTriIndex = fcount - 1;

		Surface_mesh::Face tri = mesh->getFace(randTriIndex);
		uint tri_id = tri.idx();

		// Get triangle center and normal
		Vec3d triCenter = mesh->faceCenter(tri);
		Vec3d triNor = *mesh->fn(tri);

		sp = SamplePoint(triCenter, triNor, mesh->faceArea(tri), tri_id, 1 / 3.0, 1 / 3.0);
	}

	return sp;
}

StdVector<SamplePoint> Sampler::getSamples(int numberSamples, double weight)
{
	StdVector<SamplePoint> samples(numberSamples);

	for(int i = 0; i < numberSamples; i++)
	{
		samples[i] = getSample(weight);
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

	interval = StdVector<AreaFace>(mesh->n_faces() + 1);
	interval[0] = AreaFace(0.0, mesh->getFace(0));
	int i = 0;

	// Compute mesh area in a cumulative manner
	for(int j = 0; j < (int) faceAreas.size(); j++)
	{
		interval[i+1] = AreaFace(interval[i].area + faceProbability[j], mesh->getFace(j));
		i++;
	}

	// Normalize new total mesh area
	double total = interval.back().area;
	for(int i = 0; i < (int)interval.size(); i++)
		interval[i].area /= total;

	isReady = true;
}

void Sampler::draw(const StdVector<SamplePoint> & samples)
{
	glEnable(GL_BLEND); 
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glPointSize(4.0);

	glBegin(GL_POINTS);
	uchar * rgb = new uchar[3];	
	for(int i = 0; i < (int)samples.size(); i++)
	{
		const SamplePoint * p = &samples[i];

		ColorMap::jetColorMap(rgb, p->weight, 0., 1.);	
		glColor4f(rgb[0],rgb[1],rgb[2],1);
		glNormal3dv(p->n);
		glVertex3dv(p->pos);
	}
	glEnd();
}

StdVector<SamplePoint> Sampler::getSamplesFromQSegMesh( QSegMesh* srcMesh, int numberSamples )
{
	StdVector<SamplePoint> samples;

	std::vector< Sampler > sampler;
	std::vector< double > area;
	for (int i = 0; i < srcMesh->nbSegments(); i++)
	{
		Sampler s( srcMesh->getSegment(i) );
		sampler.push_back( s );
		area.push_back( s.totalMeshArea );
	}

	double totalArea = Sum(area);
	for (int i = 0; i < srcMesh->nbSegments(); i++)
	{
		int n = (area[i] / totalArea) * numberSamples;
		std::vector< SamplePoint > samps = sampler[i].getSamples(n);
		samples.insert(samples.end(), samps.begin(), samps.end());
	}

	return samples;
}
