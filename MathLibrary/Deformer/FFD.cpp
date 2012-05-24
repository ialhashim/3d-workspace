#include "FFD.h"

FFD::FFD( QSurfaceMesh * src_mesh, FFD_FitType fit_type )
{
	this->mesh = src_mesh;

	if(!mesh) return;

	// Mesh dimensions
	width = mesh->bbmax.x() - mesh->bbmin.x();
	length = mesh->bbmax.y() - mesh->bbmin.y();
	height = mesh->bbmax.z() - mesh->bbmin.z();

	// Expand a bit so all the vertices are in (0,1)
	width += mesh->radius * 0.05;
	length += mesh->radius * 0.05;
	height += mesh->radius * 0.05;
}

void FFD::bbFit( Vec3i res )
{
	int Nx = Max(2, res.x());
	int Ny = Max(2, res.y());
	int Nz = Max(2, res.z());

	this->resolution = Vec3i(Nx, Ny, Nz);

	double dx = width / (Nx-1);
	double dy = length / (Ny-1);
	double dz = height / (Nz-1);

	Vec3d start_corner(-width/2, -length/2, -height/2);

	start_corner += mesh->center;

	// indexing
	int i = 0;

	// Nx x Ny x Nz
	pointsGridIdx = StdVector<StdVector<StdVector < int > > > 
		(Nx, StdVector< StdVector < int > >(Ny, StdVector < int >(Nz))); 

	for(int z = 0; z < Nz; z++){
		for(int y = 0; y < Ny; y++){
			for(int x = 0; x < Nx; x++){
				// Grid indexing
				pointsGridIdx[x][y][z] = i;

				// Control point position
				Vec3d p = start_corner + Vec3d(dx * x, dy * y, dz * z);

				// Add it
				points.push_back(new QControlPoint(p, i++, Vec3i(x,y,z)));
			}
		}	
	}

	// Setup local coordinate
	mP = start_corner; // this is the origin of our local frame
	mS = width * Vec3d(1,0,0);
	mT = length * Vec3d(0,1,0);
	mU = height * Vec3d(0,0,1);

	// Get copy of original mesh vertices in local coordinates
	Surface_mesh::Vertex_property<Point> mesh_points = mesh->vertex_property<Point>("v:point");
	Surface_mesh::Vertex_iterator vit, vend = mesh->vertices_end();

	meshVerticesLocal.reserve(mesh->n_vertices());

	for (vit = mesh->vertices_begin(); vit != vend; ++vit){
		meshVerticesLocal.push_back( getLocalCoordinates(mesh_points[vit]) );
	}
}

void FFD::apply()
{
	Surface_mesh::Vertex_property<Point> mesh_points = mesh->vertex_property<Point>("v:point");
	Surface_mesh::Vertex_iterator vit, vend = mesh->vertices_end();

	for (vit = mesh->vertices_begin(); vit != vend; ++vit)
	{
		int vidx = ((Surface_mesh::Vertex)vit).idx();

		mesh_points[vit] = getWorldCoordinate( deformVertexLocal(meshVerticesLocal[vidx]) );
	}
}

void FFD::fixed( Vec3i res, Vec3d location, double spacing, StdMap<int,Point> pnts )
{
	int Nx = Max(2, res.x());
	int Ny = Max(2, res.y());
	int Nz = Max(2, res.z());

	this->resolution = Vec3i(Nx, Ny, Nz);

	width = length = height = spacing;

	double dx = width / (Nx-1);
	double dy = length / (Ny-1);
	double dz = height / (Nz-1);

	Vec3d start_corner(-width/2, -length/2, -height/2);

	start_corner += location;

	// indexing
	int i = 0;

	// Nx x Ny x Nz
	pointsGridIdx = StdVector<StdVector<StdVector < int > > > 
		(Nx, StdVector< StdVector < int > >(Ny, StdVector < int >(Nz))); 

	for(int z = 0; z < Nz; z++){
		for(int y = 0; y < Ny; y++){
			for(int x = 0; x < Nx; x++){
				// Grid indexing
				pointsGridIdx[x][y][z] = i;

				// Control point position
				Vec3d p = start_corner + Vec3d(dx * x, dy * y, dz * z);

				// Add it
				points.push_back(new QControlPoint(p, i++, Vec3i(x,y,z)));
			}
		}	
	}

	// Setup local coordinate
	mP = start_corner; // this is the origin of our local frame
	mS = spacing * Vec3d(1,0,0);
	mT = spacing * Vec3d(0,1,0);
	mU = spacing * Vec3d(0,0,1);

	// Get copy of original mesh vertices in local coordinates
	for(StdMap<int,Point>::iterator it = pnts.begin(); it != pnts.end(); it++)
	{
		int vid = it->first;
		Point pos = it->second;
		fixedPointsLocal[vid] = getLocalCoordinates(pos);
	}
}

StdMap<int,Point> FFD::applyFixed()
{
	StdMap<int,Point> deformed;

	for(StdMap<int,Point>::iterator it = fixedPointsLocal.begin(); it != fixedPointsLocal.end(); it++)
	{
		deformed[it->first] = getWorldCoordinate( deformVertexLocal(it->second) );
	}

	return deformed;
}

Vec3d FFD::deformVertexLocal( const Vec3d & localPoint )
{
	Vec3d newVertex (0,0,0);

	double s = localPoint.x();
	double t = localPoint.y();
	double u = localPoint.z();

	int S = resolution.x()-1, T = resolution.y()-1, U = resolution.z()-1;

	// From Explicit definition of Bézier curves
	for (int k = 0; k <= U; k++)
	{
		int ck = Coef(U,k);
		double uk = pow(u,k) * pow(1-u, U-k);

		for (int j = 0; j <= T; j++)
		{
			int cj = Coef(T,j);
			double tj = pow(t,j) * pow(1-t, T-j);

			for (int i = 0; i <= S; i++)
			{
				int ci =  Coef(S,i);
				double si = pow(s,i) * pow(1-s, S-i);

				Vec3d controlPointPos = getLocalCoordinates(points[pointsGridIdx[i][j][k]]->pos);

				// as combination
				newVertex += ci*cj*ck*si*tj*uk * (controlPointPos);
			}
		}
	}

	return newVertex;
}

Vec3d FFD::getWorldCoordinate(const Vec3d & pLocal)
{
	return mP + pLocal.x()*mS + pLocal.y()*mT + pLocal.z()*mU;
}

Vec3d FFD::getLocalCoordinates( const Vec3d & p )
{
	Vec3d V = p;
	double s = 0, t = 0, u = 0;

	Vec3d TXU = cross(mT,mU);
	s = dot(TXU, V - mP) / dot(TXU, mS);

	Vec3d SXU = cross(mS,mU);
	t = dot(SXU, V - mP) / dot(SXU, mT);

	Vec3d SXT = cross(mS,mT);
	u = dot(SXT, V - mP) / dot(SXT, mU);

	return Vec3d(s,t,u);
}
