// From Geometric Tools, LLC

#include "ConvexHull3.h"
#include "Macros.h"

int CH_PRECISION;

ConvexHull3::ConvexHull3( std::vector<Vector3> pnts )
{
	mPnts = pnts;
	isReady = false;
	epsilon = Epsilon_LOW;

	while(!computeCH())
	{
		CH_PRECISION++;
		printf("INFO: ConvexHull precision changed to : %d \n", CH_PRECISION);

		// rest
		mPnts = pnts;
	}
}

ConvexHull3::ConvexHull3( Surface_mesh * mesh )
{
	// Get points
	Surface_mesh::Vertex_property<Point> points = mesh->vertex_property<Point>("v:point");
	Surface_mesh::Vertex_iterator vit, vend = mesh->vertices_end();

	for (vit = mesh->vertices_begin(); vit != vend; ++vit)
		mPnts.push_back(points[vit]);

	std::vector<Vector3> copyPnts = mPnts;

	isReady = false;

	// Compute CH
	while(!computeCH())
	{
		CH_PRECISION++;
		printf("INFO: ConvexHull precision changed to : %d \n", CH_PRECISION);
		mPnts = copyPnts;
	}
}


bool ConvexHull3::computeCH()
{
	std::vector<int> mExtreme;
	bool CCW = getExtremes(mExtreme);

	int i0 = mExtreme[0];
	int i1 = mExtreme[1];
	int i2 = mExtreme[2];
	int i3 = mExtreme[3];

	epsilon = (mPnts[i0]-mPnts[i1]).norm() / (Real)pow(10.0, CH_PRECISION);


	TriFace* tri0;
	TriFace* tri1;
	TriFace* tri2;
	TriFace* tri3;

	if (CCW){
		tri0 = new TriFace(i0, i1, i3);
		tri1 = new TriFace(i0, i2, i1);
		tri2 = new TriFace(i0, i3, i2);
		tri3 = new TriFace(i1, i2, i3);
		tri0->AttachTo(tri1, tri3, tri2);
		tri1->AttachTo(tri2, tri3, tri0);
		tri2->AttachTo(tri0, tri3, tri1);
		tri3->AttachTo(tri1, tri2, tri0);
	}
	else{
		tri0 = new TriFace(i0, i3, i1);
		tri1 = new TriFace(i0, i1, i2);
		tri2 = new TriFace(i0, i2, i3);
		tri3 = new TriFace(i1, i3, i2);
		tri0->AttachTo(tri2, tri3, tri1);
		tri1->AttachTo(tri0, tri3, tri2);
		tri2->AttachTo(tri1, tri3, tri0);
		tri3->AttachTo(tri0, tri2, tri1);
	}

	mHull.clear();
	mHull.insert(tri0);
	mHull.insert(tri1);
	mHull.insert(tri2);
	mHull.insert(tri3);

	for (int i=0;i<mPnts.size();i++)
	{
		if (!Update(i)){
			DeleteHull();
			return true;
		}

		if(mPnts.size() == 0)
			return false;
	}

	ExtractIndices();

	isReady = true;

	return true;
}

bool ConvexHull3::Update (int i)
{
	// Locate a TriFace visible to the input point (if possible).
    TriFace* visible = 0;
    TriFace* tri;
    std::set<TriFace*>::iterator iter = mHull.begin();
    std::set<TriFace*>::iterator end = mHull.end();
    for (/**/; iter != end; ++iter)
    {
        tri = *iter;
        if (tri->GetSign(i, mPnts, epsilon) > 0)
        {
            visible = tri;
            break;
        }
    }

    if (!visible)
    {
        // The point is inside the current hull; nothing to do.
        return true;
    }

    // Locate and remove the visible TriFaces.
    std::stack<TriFace*> visibleSet;
    std::map<int,TerminatorData> terminator;
    visibleSet.push(visible);
    visible->OnStack = true;
    int j, v0, v1;

    while (!visibleSet.empty())
    {
        tri = visibleSet.top();
        visibleSet.pop();
        tri->OnStack = false;
        for (j = 0; j < 3; ++j)
        {
            TriFace* adj = tri->Adj[j];
            if (adj)
            {
                // Detach TriFace and adjacent TriFace from each other.
                int nullIndex = tri->DetachFrom(j, adj);

				if(tri->Time != -1)
				{
					mPnts.clear();
					return true;
				}

                if (adj->GetSign(i, mPnts, epsilon) > 0)
                {
                    if (!adj->OnStack)
                    {
                        // Adjacent TriFace is visible.
                        visibleSet.push(adj);
                        adj->OnStack = true;
                    }
                }
                else
                {
                    // Adjacent TriFace is invisible.
                    v0 = tri->V[j];
                    v1 = tri->V[(j+1)%3];
                    terminator[v0] = TerminatorData(v0, v1, nullIndex, adj);
                }
            }
        }

		mHull.erase(tri);

		delete tri;
		tri = NULL;
    }

    // Insert the new edges formed by the input point and the terminator
    // between visible and invisible TriFaces.
    int size = (int)terminator.size();
    std::map<int,TerminatorData>::iterator edge = terminator.begin();

	if(edge == terminator.end())
	{
		mPnts.clear();
		return true;
	}

    v0 = edge->second.V[0];
    v1 = edge->second.V[1];
    tri = new TriFace(i, v0, v1);
    mHull.insert(tri);

    // Save information for linking first/last inserted new TriFaces.
    int saveV0 = edge->second.V[0];
    TriFace* saveTri = tri;

    // Establish adjacency links across terminator edge.
    tri->Adj[1] = edge->second.T;
    edge->second.T->Adj[edge->second.NullIndex] = tri;
    for (j = 1; j < size; ++j)
    {
        edge = terminator.find(v1);

		if(edge == terminator.end())
		{
			mPnts.clear();
			return true;
		}

        v0 = v1;
        v1 = edge->second.V[1];

        TriFace* next = new TriFace(i, v0, v1);
        mHull.insert(next);

        // Establish adjacency links across terminator edge.
        next->Adj[1] = edge->second.T;
        edge->second.T->Adj[edge->second.NullIndex] = next;

        // Establish adjacency links with previously inserted TriFace.
        next->Adj[0] = tri;
        tri->Adj[2] = next;

        tri = next;
    }

    // Establish adjacency links between first/last TriFaces.
    saveTri->Adj[0] = tri;
    tri->Adj[2] = saveTri;
    return true;
}

void ConvexHull3::ExtractIndices ()
{
	mNumSimplices = (int)mHull.size();

	std::set<TriFace*>::iterator iter = mHull.begin();
	std::set<TriFace*>::iterator end = mHull.end();
	for (; iter != end; ++iter)
	{
		TriFace* tri = *iter;
		for (int j = 0; j < 3; ++j)
		{
			mIndices.push_back(tri->V[j]);
		}
		delete tri;
		tri = NULL;
	}
	mHull.clear();
}

void ConvexHull3::DeleteHull ()
{
	std::set<TriFace*>::iterator iter = mHull.begin();
	std::set<TriFace*>::iterator end = mHull.end();
	for (/**/; iter != end; ++iter)
	{
		TriFace* tri = *iter;
		delete tri;
		tri = NULL;
	}
	mHull.clear();
}

bool ConvexHull3::getExtremes( std::vector<int> &mExtreme )
{
	bool mExtremeCCW = false;
	mExtreme.resize(4, 0);
	int nbPnts = mPnts.size();

	// Compute the axis-aligned bounding box for the input mPnts.  Keep track
	// of the indices into 'mPnts' for the current min and max.
	Vector3 mMin, mMax;
	int j, indexMin[3], indexMax[3];
	for (j = 0; j < 3; ++j)
	{
		mMin[j] = mPnts[0][j];
		mMax[j] = mMin[j];
		indexMin[j] = 0;
		indexMax[j] = 0;
	}

	int i;
	for (i = 1; i < nbPnts; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			if (mPnts[i][j] < mMin[j])
			{
				mMin[j] = mPnts[i][j];
				indexMin[j] = i;
			}
			else if (mPnts[i][j] > mMax[j])
			{
				mMax[j] = mPnts[i][j];
				indexMax[j] = i;
			}
		}
	}

	// Determine the maximum range for the bounding box.
	Real mMaxRange = mMax[0] - mMin[0];
	mExtreme[0] = indexMin[0];
	mExtreme[1] = indexMax[0];
	Real range = mMax[1] - mMin[1];
	if (range > mMaxRange)
	{
		mMaxRange = range;
		mExtreme[0] = indexMin[1];
		mExtreme[1] = indexMax[1];
	}
	range = mMax[2] - mMin[2];
	if (range > mMaxRange)
	{
		mMaxRange = range;
		mExtreme[0] = indexMin[2];
		mExtreme[1] = indexMax[2];
	}

	// The origin is either the point of minimum x-value, point of
	// minimum y-value, or point of minimum z-value.
	Vector3 mOrigin = mPnts[mExtreme[0]];


	// Test whether the point set is (nearly) a line segment.
	std::vector<Vector3> mDirection(3, Vector3());
	mDirection[0] = mPnts[mExtreme[1]] - mOrigin;
	mDirection[0].normalize();
	Real maxDistance = (Real)0;
	Real distance, Dot;
	mExtreme[2] = mExtreme[0];
	for (i = 0; i < nbPnts; ++i)
	{
		Vector3 diff = mPnts[i] - mOrigin;
		Dot = dot(mDirection[0], diff);
		Vector3 proj = diff - Dot*mDirection[0];
		distance = proj.norm();
		if (distance > maxDistance)
		{
			maxDistance = distance;
			mExtreme[2] = i;
		}
	}



	// Test whether the point set is (nearly) a planar polygon.
	mDirection[1] = mPnts[mExtreme[2]] - mOrigin;
	Dot = dot(mDirection[0], mDirection[1]);
	mDirection[1] -= Dot*mDirection[0];
	mDirection[1].normalize();
	mDirection[2] = cross(mDirection[0], mDirection[1]);
	maxDistance = (Real)0;
	Real maxSign = (Real)0;
	mExtreme[3] = mExtreme[0];
	for (i = 0; i < nbPnts; ++i)
	{
		Vector3 diff = mPnts[i] - mOrigin;
		distance = dot(mDirection[2], diff);
		if (fabs(distance) > maxDistance)
		{
			maxDistance = fabs(distance);
			mExtremeCCW = (distance > (Real)0);
			mExtreme[3] = i;
		}
	}

	return mExtremeCCW;
}

void ConvexHull3::draw()
{
	if (!isReady) return;

	std::vector<std::vector<Vector3>> tris;
	for (int i=0;i<3*mNumSimplices;)
	{
		std::vector<Vector3> tri;

		tri.push_back(mPnts[mIndices[i++]]);
		tri.push_back(mPnts[mIndices[i++]]);
		tri.push_back(mPnts[mIndices[i++]]);

		tris.push_back(tri);
	}

	SimpleDraw::DrawTriangles(tris);
}

// ConvexHull3::TriFace
ConvexHull3::TriFace::TriFace (int v0, int v1, int v2)
: Sign(0),Time(-1),OnStack(false)
{
	V[0] = v0;
	V[1] = v1;
	V[2] = v2;
	Adj[0] = nullptr;
	Adj[1] = nullptr;
	Adj[2] = nullptr;
}

int ConvexHull3::TriFace::GetSign( int id , std::vector<Vector3> &pnts, Real epsilon )
{
	// check indices
	if(V[0] < 0 || V[0] >= pnts.size()) return 1;
	if(V[1] < 0 || V[1] >= pnts.size()) return 1;
	if(V[2] < 0 || V[2] >= pnts.size()) return 1;

	Vector3 ab = pnts[V[1]] - pnts[V[0]];
	ab.normalize();
	Vector3 ac = pnts[V[2]] - pnts[V[0]];
	ac.normalize();
	Vector3 cross_prod = cross(ab,ac);
	cross_prod.normalize();
	Vector3 av = pnts[id]	- pnts[V[0]];
	return dot(cross_prod, av) > epsilon;
}



void ConvexHull3::TriFace::AttachTo (TriFace* adj0, TriFace* adj1,
	TriFace* adj2)
{
	// assert:  The input adjacent TriFaces are correctly ordered.
	Adj[0] = adj0;
	Adj[1] = adj1;
	Adj[2] = adj2;
}

int ConvexHull3::TriFace::DetachFrom (int adjIndex, TriFace* adj)
{
	Adj[adjIndex] = nullptr;
	for (int i = 0; i < 3; ++i)
	{
		if (adj->Adj[i] == this)
		{
			adj->Adj[i] = nullptr;
			return i;
		}
	}
	return -1;
}

// ConvexHull3::TerminatorData
ConvexHull3::TerminatorData::TerminatorData (int v0, int v1, int nullIndex, TriFace* tri) 
	: NullIndex(nullIndex), T(tri)
{
	V[0] = v0;
	V[1] = v1;
}
