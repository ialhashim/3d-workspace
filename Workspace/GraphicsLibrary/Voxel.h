#pragma once

struct Voxel{ 
	int x, y, z;
	int flag;

	Voxel(){ x = y = z = flag = 0; }
	Voxel(int X, int Y, int Z) : x(X), y(Y), z(Z){ flag = 1; } 

	// Operators
	operator const Vec() const{ return Vec(x,y,z); }

	Voxel & operator+= (const Voxel & other){
		x += other.x;	y += other.y;	z += other.z;
		return *this;
	}

	Voxel operator+(const Voxel & other) const{
		return Voxel(*this) += other;
	}

	// Useful for bounds
	inline void toMax(const Voxel & v){ x = Max(x, v.x); y = Max(y, v.y); z = Max(z, v.z); }
	inline void toMin(const Voxel & v){ x = Min(x, v.x); y = Min(y, v.y); z = Min(z, v.z); }
};
