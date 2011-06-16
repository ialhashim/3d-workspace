#pragma once

#include "Matrixf.h"

#define M_PI 3.14159265358979323846

// Noise code based on Ken Perlin's Improved Noise http://mrl.nyu.edu/~perlin/noise/
static int p[] = { 151,160,137,91,90,15,
	131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
	190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
	88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
	77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
	102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
	135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
	5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
	223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
	129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
	251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
	49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
	138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180,
	151,160,137,91,90,15,
	131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
	190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
	88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
	77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
	102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
	135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
	5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
	223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
	129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
	251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
	49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
	138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180
};

class NoiseGen
{
private:

	static double noise(double x, double y, double z){
		int X = (int)floor(x) & 255, Y = (int)floor(y) & 255, Z = (int)floor(z) & 255;

		x -= floor(x); y -= floor(y); z -= floor(z);

		double u = fade(x), v = fade(y), w = fade(z);

		int A = p[X]+Y, AA = p[A]+Z, AB = p[A+1]+Z, B = p[X+1]+Y, BA = p[B]+Z, BB = p[B+1]+Z; 

		return lerp(w, lerp(v, lerp(u, grad(p[AA  ], x  , y  , z   ),
			grad(p[BA  ], x-1, y  , z   )),
			lerp(u, grad(p[AB  ], x  , y-1, z   ), 
			grad(p[BB  ], x-1, y-1, z   ))),
			lerp(v, lerp(u, grad(p[AA+1], x  , y  , z-1 ),
			grad(p[BA+1], x-1, y  , z-1 )), 
			lerp(u, grad(p[AB+1], x  , y-1, z-1 ),
			grad(p[BB+1], x-1, y-1, z-1 ))));
	}

	static double fade(double t) { return t * t * t * (t * (t * 6 - 15) + 10); }
	static double lerp(double t, double a, double b) { return a + t * (b - a); }
	static double grad(int hash, double x, double y, double z) {
		int h = hash & 15;                     
		double u = h<8 ? x : y, v = h<4 ? y : h==12||h==14 ? x : z;
		return ((h&1) == 0 ? u : -u) + ((h&2) == 0 ? v : -v);
	}

	static float period(float t)
	{
		return cos(M_PI * t);
	}

public:

	static MatrixXf noise(int width, int height)
	{
		MatrixXf img = MatrixXf::Random(height, width);

		// Normalize
		img.array() = (img.array() + img.array().minCoeff());
		img.array() /= img.array().maxCoeff();

		return img;
	}

	static MatrixXf smoothNoise(int width, int height, float scale = 1.0f)
	{
		MatrixXf img = MatrixXf::Zero(height, width);

		float MinDimension = Min(width, height) * scale;

		for(int y = 0; y < height; y++)
		{
			for(int x = 0; x < width; x++)
			{
				double u = (float)x / MinDimension;
				double v = (float)y / MinDimension;

				v = abs(v - 0.5 * (1/scale));

				img(y,x) = (float)noise(u,v,1);
			}
		}

		return img;
	}

	static MatrixXf wavePattern(int width, int height, int cycleCount = 4, bool vertical = false)
	{
		MatrixXf img = MatrixXf::Zero(height, width);
		VectorXf v = VectorXf::Zero(width);

		int len = width / cycleCount;

		for(int c = 0; c < cycleCount + 1; c++)
		{
			for(int i = 0; i < len; i++)
			{
				float t = (c+1 % 2) - (float)i / len;

				int index = i+(c*len);

				if(index < width)
					v(index) = period(t);
				else
					break;
			}
		}

		img = v.rowwise().replicate(height);

		img.array() = img.array() * 0.5f;

		if(vertical)
			return img.transpose();
		else
			return img;
	}
};
