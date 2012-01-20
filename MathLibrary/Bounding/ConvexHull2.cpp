// Geometric Tools, LLC

#include "ConvexHull2.h"


ConvexHull2::ConvexHull2(std::vector<Vector2> &pnts)
	:mVertices(pnts)
{
	std::vector<int> mExtreme;
	bool CCW = getExtremes(mExtreme);
    int i0 = mExtreme[0];
    int i1 = mExtreme[1];
    int i2 = mExtreme[2];


    int i;
	int mNumVertices = mVertices.size();
    Edge* edge0;
    Edge* edge1;
    Edge* edge2;

    if (CCW)
    {
        edge0 = new Edge(i0, i1);
        edge1 = new Edge(i1, i2);
        edge2 = new Edge(i2, i0);
    }
    else
    {
        edge0 = new Edge(i0, i2);
        edge1 = new Edge(i2, i1);
        edge2 = new Edge(i1, i0);
    }

    edge0->Insert(edge2, edge1);
    edge1->Insert(edge0, edge2);
    edge2->Insert(edge1, edge0);

    Edge* hull = edge0;
    for (i = 0; i < mNumVertices; ++i)
    {
        if (!Update(hull, i))
        {
            hull->DeleteAll();
            return;
        }
    }

    hull->GetIndices(mNumSimplices, mIndices);
    hull->DeleteAll();
}

bool ConvexHull2::getExtremes( std::vector<int> &mExtreme )
{
	bool mExtremeCCW = false;
	mExtreme.resize(3, 0);
	int nummVertices = mVertices.size();

	// Compute the axis-aligned bounding box for the input mVertices.  Keep track
	// of the indices into 'mVertices' for the current min and max.
	int j, indexMin[2], indexMax[2];
	Vector2 mMin, mMax;
	for (j = 0; j < 2; ++j)
	{
		mMin[j] = mVertices[0][j];
		mMax[j] = mMin[j];
		indexMin[j] = 0;
		indexMax[j] = 0;
	}

	int i;
	for (i = 1; i < nummVertices; ++i)
	{
		for (j = 0; j < 2; ++j)
		{
			if (mVertices[i][j] < mMin[j])
			{
				mMin[j] = mVertices[i][j];
				indexMin[j] = i;
			}
			else if (mVertices[i][j] > mMax[j])
			{
				mMax[j] = mVertices[i][j];
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

	// The origin is either the point of minimum x-value or point of
	// minimum y-value.
	Vector2 mOrigin = mVertices[mExtreme[0]];


	// Test whether the point set is (nearly) a line segment.
	std::vector<Vector2> mDirection(2, Vector2());
	mDirection[0] = mVertices[mExtreme[1]] - mOrigin;
	mDirection[0].normalize();
	mDirection[1] = Vector2(-mDirection[0].y(), mDirection[0].x());
	Real maxDistance = (Real)0;
	mExtreme[2] = mExtreme[0];
	for (i = 0; i < nummVertices; ++i)
	{
		Vector2 diff = mVertices[i] - mOrigin;
		Real distance = dot(mDirection[1], diff);

		if (fabs(distance) > maxDistance)
		{
			maxDistance = fabs(distance);
			mExtremeCCW = (distance > (Real)0);
			mExtreme[2] = i;
		}
	}


	return mExtremeCCW;
}

bool ConvexHull2::Update (Edge*& hull, int i)
{
    // Locate an edge visible to the input point (if possible).
    Edge* visible = 0;
    Edge* current = hull;
    do
    {
        if (current->GetSign(i, mVertices) > 0)
        {
            visible = current;
            break;
        }

        current = current->E[1];
    }
    while (current != hull);

    if (!visible)
    {
        // The point is inside the current hull; nothing to do.
        return true;
    }

    // Remove the visible edges.
    Edge* adj0 = visible->E[0];
    if (!adj0)
    {
        return false;
    }

    Edge* adj1 = visible->E[1];
    if (!adj1)
    {
        return false;
    }

    visible->DeleteSelf();

    while (adj0->GetSign(i, mVertices) > 0)
    {
        hull = adj0;
        adj0 = adj0->E[0];
        if (!adj0)
        {
            return false;
        }

        adj0->E[1]->DeleteSelf();
    }

    while (adj1->GetSign(i,mVertices) > 0)
    {
        hull = adj1;
        adj1 = adj1->E[1];
        if (!adj1)
        {
            return false;
        }

        adj1->E[0]->DeleteSelf();
    }

    // Insert the new edges formed by the input point and the end mVertices of
    // the polyline of invisible edges.
    Edge* edge0 = new Edge(adj0->V[1], i);
    Edge* edge1 = new Edge(i, adj1->V[0]);
    edge0->Insert(adj0, edge1);
    edge1->Insert(edge0, adj1);
    hull = edge0;

    return true;
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// ConvexHull2::Edge
//----------------------------------------------------------------------------

ConvexHull2::Edge::Edge (int v0, int v1)
    :Sign(0),  Time(-1)
{
    V[0] = v0;
    V[1] = v1;
    E[0] = 0;
    E[1] = 0;
}
//----------------------------------------------------------------------------

int ConvexHull2::Edge::GetSign( int i, std::vector<Vector2> &pnts )
{
	Vector2 u = pnts[V[1]]  - pnts[V[0]];
	Vector2 v = Vector2(-u.y(), u.x());
	Vector2 p = pnts[i] - pnts[V[0]];
	
	return dot(v,p)<0;
}
//----------------------------------------------------------------------------

void ConvexHull2::Edge::Insert (Edge* adj0, Edge* adj1)
{
    adj0->E[1] = this;
    adj1->E[0] = this;
    E[0] = adj0;
    E[1] = adj1;
}
//----------------------------------------------------------------------------

void ConvexHull2::Edge::DeleteSelf ()
{
    if (E[0])
    {
        E[0]->E[1] = 0;
    }

    if (E[1])
    {
        E[1]->E[0] = 0;
    }

    Edge* tmpThis = this;
    delete tmpThis;
}
//----------------------------------------------------------------------------

void ConvexHull2::Edge::DeleteAll ()
{
    Edge* adj = E[1];
    while (adj && adj != this)
    {
        Edge* save = adj->E[1];
        delete adj;
        adj = save;
    }

    Edge* tmpThis = this;
    delete tmpThis;
}
//----------------------------------------------------------------------------
void ConvexHull2::Edge::GetIndices( int& numIndices, std::vector<int> &indices )
{
    // Count the number of edge vertices and allocate the index array.
    numIndices = 0;
    Edge* current = this;
    do
    {
        ++numIndices;
        current = current->E[1];
    }
    while (current != this);
    indices.resize(numIndices);

    // Fill the index array.
    numIndices = 0;
    current = this;
    do
    {
        indices[numIndices] = current->V[0];
        ++numIndices;
        current = current->E[1];
    }
    while (current != this);
}



