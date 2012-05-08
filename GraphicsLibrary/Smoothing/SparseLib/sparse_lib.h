#undef Vector

#define COMPLEX std::complex<double>

#include "compcol_double.h"
#include "mvblasd.h"

// Pre conditioners
#include "diagpre_double.h"
#include "icpre_double.h"
#include "ilupre_double.h"

// Solvers
#include "bicg.h"
#include "bicgstab.h"
#include "cg.h"
#include "cgs.h"

#include "cg.hpp"
#include "mvmd.hpp"
#include "mvvd.hpp"
#include "qsort_double.hpp"
#include "spmm.hpp"
#include "spsm.hpp"

#include "compcol_double.hpp"
#include "mvblasd.hpp"

// Preconditioners
#include "diagpre_double.hpp"
#include "icpre_double.hpp"
#include "ilupre_double.hpp"

#include "mvvi.h"
#include "mvvi.hpp"