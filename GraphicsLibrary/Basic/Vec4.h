struct Vec4{
	double v[4];

	const double &operator [] (int i) const
	{ return v[i]; }
	double &operator [] (int i)
	{ return v[i]; }

	Vec4() {v[0] = v[1] = v[2] = v[3];}
	Vec4(const double from[4])
	{
		for(int i = 0; i < 4; i++)
			v[i] = from[i];
	}

	Vec4 &operator += (const Vec4 &x)
	{
		for (int i = 0; i < 4; i++)
		{
			#pragma omp atomic
			v[i] += x[i];
		}
		return *this;
	}
};

// Vec/scalar operators
static inline const Vec4 operator * (const double &x, const Vec4 &v)
{
	Vec4 result;
	for (int i = 0; i < 4; i++)
		result[i] = x * v[i];
	return result;
}
