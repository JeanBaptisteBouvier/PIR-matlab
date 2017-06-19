/*=========================================================================
 * ode78_bcp_rel.h header mex file to generate
 * a Runge-Kutta 7/8 integrator of the CR3BP vector field.
 * 
 * author:  BLB
 * year:    2015
 * version: 1.0
 * author:  JBB
 * year:    2017
 * version: 2.0
 *=======================================================================*/

//-------------------------------------------------------------------------
// Headers
//-------------------------------------------------------------------------
//Mex
#include "mex.h"

//Custom
#include "ode78_rel.h"


//-------------------------------------------------------------------------
// The gateway function.
// The input must be, in that order:
// 1. [t0, tf] the time span
// 2. double y0[12 or 48], the initial state
// 3. double mu, the cr3bp mass ratio
// 4. double theta0, the initial phase of the Sun
// 5. double ms, the mass of the Sun
// 6. double as, the semi-major axis of the Sun
// 7. double omS, the mean motion of the Sun
//-------------------------------------------------------------------------
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]);
