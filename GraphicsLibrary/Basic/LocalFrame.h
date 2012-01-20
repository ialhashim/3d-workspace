#pragma once

#include "Plane.h"

class LocalFrame
{

private:

public:
        Vec n, up, b; // should be private ?

        LocalFrame();
        LocalFrame(const LocalFrame& from);
        LocalFrame& operator= (const LocalFrame& from);
        LocalFrame(const Vec& normal, const Vec& VecUp);
        LocalFrame(const LocalFrame& prev, Vec newNormal);

        static Vector<LocalFrame> alongTangent( const Vector<Vec>& tangent );
        static Vector<LocalFrame> alongTangent( const Vector<Vec>& tangent, const Vec& firstUp );

        void draw(const Vec& center);

        static qglviewer::Quaternion rotation(const LocalFrame& A, const LocalFrame& B);
};

