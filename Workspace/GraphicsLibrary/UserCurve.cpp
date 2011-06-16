#include "UserCurve.h"

UserCurve::UserCurve()
{
	isReady = false;
	isSimplified = false;
}

Vec UserCurve::intersectionRayPlane(const Vec & planeNormal, const Vec & planeOrigin, const Vec & rayOrigin, const Vec & rayDirection)
{
	float dot = planeNormal * rayDirection;

	if(dot >= 0)
		return Vec();
	else
	{
		float distance = ( planeNormal * (planeOrigin - rayOrigin))  / (planeNormal * rayDirection);
		return rayOrigin + (rayDirection * distance);
	}
}

void UserCurve::addPoint(const Vec & origin, const Vec & direction)
{
	if(!(plane_normal.x == 0 && plane_normal.y == 0 && plane_normal.z == 0))
	{
		// find intersection point 'p'
		Vec p = intersectionRayPlane(plane_normal.unit(), plane_origin, origin, direction.unit());

		if(p.norm() > 0)
		{
			points.push_back(p);

			if(spline.GetNumPoints() > 4)
			{
				float d_last_point = (p - spline.GetNthPoint(spline.GetNumPoints() - 1)).norm();

				if(d_last_point > 0.005f)
				{
					spline.AddSplinePoint(p);
					isReady = true;
				}
			}
			else
			{
				spline.AddSplinePoint(p);
			}
		}
	}
}

void UserCurve::setPlane(const Vec & pointOnPlane, const Vec & normal)
{
	plane_origin = pointOnPlane;
	plane_normal = normal.unit();

	// Reset curve
	points.clear();
	spline.Clear();

	isSimplified = false;
	isReady = false;

	isVisible = true;
}

Plane UserCurve::getPlane()
{
	return Plane(plane_normal, plane_origin);
}

void UserCurve::moveHead(const Vec & toHead)
{
	Vec delta = toHead - spline.GetNthPoint(0);

	spline.Translate(delta);
}

void UserCurve::simplify()
{
	if(spline.GetNumPoints() > 8)
	{
		spline.Crop(4, 1);
		spline.Simplify(1);
		spline.Smooth(1);
		spline.Subdivide(3);
	}

	isSimplified = true;
}

void UserCurve::simplify2()
{
	spline.Smooth(2);
	spline.Subdivide(3);
	
	isSimplified = true;
}

Vector<Vec> UserCurve::getPath(float stepLength)
{
	return spline.GetUniformPath(stepLength);
}

Spline * UserCurve::getSpline()
{
	return &spline;
}

void UserCurve::draw()
{
	if(isVisible == true)
	{
		// create plane indicator
		Vec p1, p2, p3, p4;

		Vec u = plane_normal.orthogonalVec().unit();
		Vec v = (u ^ plane_normal).unit();

		float width = 0.025f;

		p1 = u * width + -v * width;
		p2 = u * width + v * width;
		p4 = -u * width + -v * width;
		p3 = -u * width + v * width;

		glClear(GL_DEPTH_BUFFER_BIT);
		glDisable(GL_LIGHTING);

		glLineWidth(2);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

		glDisable (GL_BLEND);

		// blueish color
		glColor3f(0,0.5,1);

		// Draw plane
		glBegin(GL_QUADS);
			glVertex3fv(p1 + plane_origin);
			glVertex3fv(p2 + plane_origin);
			glVertex3fv(p3 + plane_origin);
			glVertex3fv(p4 + plane_origin);
		glEnd();

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		//glPointSize(5);
		//glColor4f(1,1,0,0.25);
		// Draw points
		/*glBegin(GL_POINTS);
		for(int i = 0; i < points.size(); i++)
			glVertex3fv(points[i]);
		glEnd();*/

		// Setup spline render options

		glLineWidth(5);
		glColor3f(1,1,0);

		glEnable(GL_BLEND); 
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		Vector<Vec> pointsDraw, smoothed;
		for(int i = 0; i < spline.GetNumPoints(); i++)
			pointsDraw.push_back(spline.GetNthPoint(i));

		// smooth points to draw
		smoothed = pointsDraw;
                for(int i = 1; i < (int)pointsDraw.size() - 1; i++)
			smoothed[i] = (pointsDraw[i-1] + pointsDraw[i+1]) / 2.0f;
		pointsDraw = smoothed;
		
		if(spline.GetNumPoints() > 3)
		{
			// Draw Spline
			glBegin(GL_LINE_STRIP);
                        for(int i = 0; i < (int)pointsDraw.size(); i++)
				glVertex3fv(pointsDraw[i]);
			glEnd();

			/*glPointSize(8);
			glColor3f(0,0,1);
			glBegin(GL_POINTS);
			for(int i = 0; i < spline.GetNumPoints(); i++)
				glVertex3fv(spline.GetNthPoint(i));
			glEnd();*/
		}

	}
}
