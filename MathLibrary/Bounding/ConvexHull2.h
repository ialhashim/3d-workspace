// Geometric Tools, LLC

#include "GraphicsLibrary/Mesh/SurfaceMesh/Vector.h"
#include <vector>


typedef double								Real;
typedef Vector<Real, 2>						Vector2;

class  ConvexHull2 
{
public:
    ConvexHull2 (std::vector<Vector2> &pnts);
	int getNumSimplices(){return mNumSimplices;}
	std::vector<int>& getIndices(){return mIndices;}

private:
	class Edge
	{
	public:
		Edge (int v0, int v1);

		int GetSign (int i, std::vector<Vector2> &pnts);
		void Insert (Edge* adj0, Edge* adj1);
		void DeleteSelf ();
		void DeleteAll ();
		void GetIndices (int& numIndices, std::vector<int> &indices);

		int V[2];
		Edge* E[2];
		int Sign;
		int Time;
	};

private:
   std::vector<Vector2> mVertices;		// The input points.
   int					mNumSimplices;
   std::vector<int>		mIndices;

private:
   bool getExtremes( std::vector<int> &mExtremes);
   bool Update (Edge*& hull, int i);
};
