#include "LocalFrame.h"
#include "SimpleDraw.h"

LocalFrame::LocalFrame()
{
}

LocalFrame::LocalFrame( const Vec& normal, const Vec& VecUp )
{
	this->n = normal;
	this->up = VecUp;
	this->b = n ^ up;
}

LocalFrame::LocalFrame( const LocalFrame& prev, Vec newNormal )
{
        qglviewer::Quaternion q(prev.n, newNormal);

	this->n = newNormal;
	this->up = q.rotate(prev.up);
	this->b = q.rotate(prev.b);
}

LocalFrame::LocalFrame( const LocalFrame& from )
{
	this->n = from.n;
	this->up = from.up;
	this->b = from.b;
}

LocalFrame& LocalFrame::operator= (const LocalFrame& from)
{
        this->n = from.n;
        this->up = from.up;
        this->b = from.b;
        return *this;
}

void LocalFrame::draw(const Vec& center)
{
	SimpleDraw::IdentifyArrow(center, center + (up * 0.01), 4);
	SimpleDraw::IdentifyArrow(center, center + (b * 0.01));
}

Vector<LocalFrame> LocalFrame::alongTangent( const Vector<Vec>& tangent )
{
        Vector<LocalFrame> result;

	// First frame
        result.push_back(LocalFrame (tangent.front(), tangent.front().orthogonalVec()));

	// Remaining
        for(unsigned int i = 1; i < tangent.size(); i++)
	{
            Vec newNormal = tangent[i];
            LocalFrame frame;
            LocalFrame prev = result.back();

            qglviewer::Quaternion q(prev.n, newNormal);

            frame.n = newNormal;
            frame.up = q.rotate(prev.up);
            frame.b = q.rotate(prev.b);

            result.push_back(frame);
        }

        return result;
}

Vector<LocalFrame> LocalFrame::alongTangent( const Vector<Vec>& tangent, const Vec& firstUp )
{
        Vector<LocalFrame> result;

	// First frame
	result.push_back(LocalFrame(tangent.front(), firstUp));

        for(int i = 1; i < (int)tangent.size(); i++)
		result.push_back( LocalFrame(result.back(), tangent[i]) );

	return result;
}

qglviewer::Quaternion LocalFrame::rotation( const LocalFrame& A, const LocalFrame& B )
{
	/*Vector<Vec> axis1, axis2;

	axis1.push_back(A.n); axis1.push_back(A.up); axis1.push_back(A.b);
	axis2.push_back(B.n); axis2.push_back(B.up); axis2.push_back(B.b);

	Transformation t = Transform3D::findFrom(axis1, axis2);*/

	// Warning
	return qglviewer::Quaternion();
}
