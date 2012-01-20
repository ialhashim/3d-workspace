// Geometric Tools, LLC
#include "MinOBB2.h"

MinOBB2::MinOBB2 (std::vector<Vector2> &points)
{
    // Get the convex hull of the points.
    std::vector<Vector2> hullPoints;

	ConvexHull2 hull(points);
	int hullNumSimplices = hull.getNumSimplices();
	std::vector<int> hullIndices = hull.getIndices();

	
	int numPoints = hullNumSimplices;
	hullPoints.resize(numPoints);
	for (int i = 0; i < numPoints; ++i)
	{
		hullPoints[i] = points[hullIndices[i]];
	}

    
    // The input points are V[0] through V[N-1] and are assumed to be the
    // vertices of a convex polygon that are counterclockwise ordered.  The
    // input points must not contain three consecutive collinear points.

    // Unit-length edge directions of convex polygon.  These could be
    // precomputed and passed to this routine if the application requires it.
    int numPointsM1 = numPoints -1;
    std::vector<Vector2> edges(numPoints, Vector2());
    std::vector<bool> visited(numPoints, false);
    int i;
    for (i = 0; i < numPointsM1; ++i)
    {
        edges[i] = hullPoints[i + 1] - hullPoints[i];
        edges[i].normalize();
        visited[i] = false;
    }
    edges[numPointsM1] = hullPoints[0] - hullPoints[numPointsM1];
    edges[numPointsM1].normalize();
    visited[numPointsM1] = false;

    // Find the smallest axis-aligned box containing the points.  Keep track
    // of the extremum indices, L (left), R (right), B (bottom), and T (top)
    // so that the following constraints are met:
    //   V[L].X() <= V[i].X() for all i and V[(L+1)%N].X() > V[L].X()
    //   V[R].X() >= V[i].X() for all i and V[(R+1)%N].X() < V[R].X()
    //   V[B].Y() <= V[i].Y() for all i and V[(B+1)%N].Y() > V[B].Y()
    //   V[T].Y() >= V[i].Y() for all i and V[(T+1)%N].Y() < V[T].Y()
    Real xmin = hullPoints[0].x(), xmax = xmin;
    Real ymin = hullPoints[0].y(), ymax = ymin;
    int LIndex = 0, RIndex = 0, BIndex = 0, TIndex = 0;
    for (i = 1; i < numPoints; ++i)
    {
        if (hullPoints[i].x() <= xmin)
        {
            xmin = hullPoints[i].x();
            LIndex = i;
        }
        if (hullPoints[i].x() >= xmax)
        {
            xmax = hullPoints[i].x();
            RIndex = i;
        }

        if (hullPoints[i].y() <= ymin)
        {
            ymin = hullPoints[i].y();
            BIndex = i;
        }
        if (hullPoints[i].y() >= ymax)
        {
            ymax = hullPoints[i].y();
            TIndex = i;
        }
    }

    // Apply wrap-around tests to ensure the constraints mentioned above are
    // satisfied.
    if (LIndex == numPointsM1)
    {
        if (hullPoints[0].x() <= xmin)
        {
            xmin = hullPoints[0].x();
            LIndex = 0;
        }
    }

    if (RIndex == numPointsM1)
    {
        if (hullPoints[0].x() >= xmax)
        {
            xmax = hullPoints[0].x();
            RIndex = 0;
        }
    }

    if (BIndex == numPointsM1)
    {
        if (hullPoints[0].y() <= ymin)
        {
            ymin = hullPoints[0].y();
            BIndex = 0;
        }
    }

    if (TIndex == numPointsM1)
    {
        if (hullPoints[0].y() >= ymax)
        {
            ymax = hullPoints[0].y();
            TIndex = 0;
        }
    }

    // The dimensions of the axis-aligned box.  The extents store width and
    // height for now.
    mMinBox.Center.x() = ((Real)0.5)*(xmin + xmax);
    mMinBox.Center.y() = ((Real)0.5)*(ymin + ymax);
    mMinBox.Axis[0] = Vector2(1.,0.);
    mMinBox.Axis[1] = Vector2(0.,1.);
    mMinBox.Extent[0] = ((Real)0.5)*(xmax - xmin);
    mMinBox.Extent[1] = ((Real)0.5)*(ymax - ymin);
    Real minAreaDiv4 = mMinBox.Extent[0]*mMinBox.Extent[1];

    // The rotating calipers algorithm.
    Vector2 U = Vector2(1.,0.);
    Vector2 V = Vector2(0.,1.);

    bool done = false;
    while (!done)
    {
        // Determine the edge that forms the smallest angle with the current
        // box edges.
        int flag = F_NONE;
        Real maxDot = (Real)0;

        Real Dot = dot(U, edges[BIndex]);
        if (Dot > maxDot)
        {
            maxDot = Dot;
            flag = F_BOTTOM;
        }

        Dot = dot(V, edges[RIndex]);
        if (Dot > maxDot)
        {
            maxDot = Dot;
            flag = F_RIGHT;
        }

        Dot = -dot(U, edges[TIndex]);
        if (Dot > maxDot)
        {
            maxDot = Dot;
            flag = F_TOP;
        }

        Dot = -dot(V, edges[LIndex]);
        if (Dot > maxDot)
        {
            maxDot = Dot;
            flag = F_LEFT;
        }

        switch (flag)
        {
        case F_BOTTOM:
            if (visited[BIndex])
            {
                done = true;
            }
            else
            {
                // Compute box axes with E[B] as an edge.
                U = edges[BIndex];
                V = Vector2(-U.y(), U.x());
                UpdateBox(hullPoints[LIndex], hullPoints[RIndex],
                    hullPoints[BIndex], hullPoints[TIndex], U, V,
                    minAreaDiv4);

                // Mark edge visited and rotate the calipers.
                visited[BIndex] = true;
                if (++BIndex == numPoints)
                {
                    BIndex = 0;
                }
            }
            break;
        case F_RIGHT:
            if (visited[RIndex])
            {
                done = true;
            }
            else
            {
                // Compute box axes with E[R] as an edge.
                V = edges[RIndex];
                U = Vector2(-V.y(), V.x());
                UpdateBox(hullPoints[LIndex], hullPoints[RIndex],
                    hullPoints[BIndex], hullPoints[TIndex], U, V,
                    minAreaDiv4);

                // Mark edge visited and rotate the calipers.
                visited[RIndex] = true;
                if (++RIndex == numPoints)
                {
                    RIndex = 0;
                }
            }
            break;
        case F_TOP:
            if (visited[TIndex])
            {
                done = true;
            }
            else
            {
                // Compute box axes with E[T] as an edge.
                U = -edges[TIndex];
                V = Vector2(-U.y(), U.x());
                UpdateBox(hullPoints[LIndex], hullPoints[RIndex],
                    hullPoints[BIndex], hullPoints[TIndex], U, V,
                    minAreaDiv4);

                // Mark edge visited and rotate the calipers.
                visited[TIndex] = true;
                if (++TIndex == numPoints)
                {
                    TIndex = 0;
                }
            }
            break;
        case F_LEFT:
            if (visited[LIndex])
            {
                done = true;
            }
            else
            {
                // Compute box axes with E[L] as an edge.
                V = -edges[LIndex];
                U = Vector2(-V.y(), V.x());
                UpdateBox(hullPoints[LIndex], hullPoints[RIndex],
                    hullPoints[BIndex], hullPoints[TIndex], U, V,
                    minAreaDiv4);

                // Mark edge visited and rotate the calipers.
                visited[LIndex] = true;
                if (++LIndex == numPoints)
                {
                    LIndex = 0;
                }
            }
            break;
        case F_NONE:
            // The polygon is a rectangle.
            done = true;
            break;
        }
    }

}


void MinOBB2::UpdateBox (const Vector2& LPoint,
    const Vector2& RPoint, const Vector2& BPoint,
    const Vector2& TPoint, const Vector2& U,
    const Vector2& V, Real& minAreaDiv4)
{
    Vector2 RLDiff = RPoint - LPoint;
    Vector2 TBDiff = TPoint - BPoint;
    Real extent0 = ((Real)0.5)*dot(U, RLDiff);
    Real extent1 = ((Real)0.5)*dot(V, TBDiff);
    Real areaDiv4 = extent0*extent1;
    if (areaDiv4 < minAreaDiv4)
    {
        minAreaDiv4 = areaDiv4;
        mMinBox.Axis[0] = U;
        mMinBox.Axis[1] = V;
        mMinBox.Extent[0] = fabs(extent0);
        mMinBox.Extent[1] = fabs(extent1);
        Vector2 LBDiff = LPoint - BPoint;
        mMinBox.Center = LPoint + U*extent0 + V*(extent1 - dot(V, LBDiff));
    }
}

MinOBB2::Box2 MinOBB2::getBox2()
{
	return mMinBox;
}

