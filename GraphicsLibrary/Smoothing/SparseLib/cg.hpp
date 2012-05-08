#pragma once
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                   MV++ Numerical Matrix/Vector C++ Library                */
/*                             MV++ Version 1.5                              */
/*                                                                           */
/*                                  R. Pozo                                  */
/*               National Institute of Standards and Technology              */
/*                                                                           */
/*                                  NOTICE                                   */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that this permission notice appear in all copies and             */
/* supporting documentation.                                                 */
/*                                                                           */
/* Neither the Institution (National Institute of Standards and Technology)  */
/* nor the author makes any representations about the suitability of this    */
/* software for any purpose.  This software is provided ``as is''without     */
/* expressed or implied warranty.                                            */
/*                                                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "cg.h"

#include "time.h"

namespace Solver
{
	template < class MMatrix, class MVector, class Preconditioner, class Real >
	int CG(const MMatrix &A, MVector &x, const MVector &b, const Preconditioner &M, int &max_iter, Real &tol)
	{
		Real resid;
		MVector p, z, q;
		MVector alpha(1), beta(1), rho(1), rho_1(1);

		MVector r = b - A*x;

		Real normb = norm(b);

		if (normb == 0.0) normb = 1;

		if ((resid = norm(r) / normb) <= tol){
			tol = resid;
			max_iter = 0;
			return 0;
		}

		for (int i = 1; i <= max_iter; i++) 
		{
			// Assign Z
			z = M.solve(r);

			rho.p_[0] = dot(r, z);

			// Assign P
			if (i == 1)
				p = z;
			else {
				beta.p_[0] = rho.p_[0] / rho_1.p_[0];
				p = z + beta.p_[0] * p;
			}

			// Assign Q
			q = A*p;

			alpha.p_[0] = rho.p_[0] / dot(p, q);

			// Change X and R
			x += alpha.p_[0] * p;
			r -= alpha.p_[0] * q;

			// Check tol
			if ((resid = norm(r) / normb) <= tol)
			{
				tol = resid;
				max_iter = i;
				return 0;     
			}

			rho_1.p_[0] = rho.p_[0];
		}

		tol = resid;
		return 1;
	}
}

/*#include "compcol_double.h"
#include "mvblasd.h"

// Preconditioners
#include "diagpre_double.h"
#include "icpre_double.h"
#include "ilupre_double.h"

template int  Solver::CG<class CompCol_Mat_double,class 
MV_Vector_double,class DiagPreconditioner_double,double>
(class CompCol_Mat_double const &,class MV_Vector_double &,class 
MV_Vector_double const &,class DiagPreconditioner_double const &,
int &,double &);
*/