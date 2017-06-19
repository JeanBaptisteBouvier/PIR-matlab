/*=========================================================================
 * ode78_bcp_event_rel.h header mex file to generate
 * a Runge-Kutta 7/8 integrator of the CR3BP vector field.
 * 
 * author:  BLB
 * year:    2015
 * version: 1.0
 *=======================================================================*/

//-------------------------------------------------------------------------
// Headers
//-------------------------------------------------------------------------
//Mex
#include "mex.h"
#include <string.h>
//Custom
#include "ode78_rel.h"
#include "custom_ode.h"
#include "custom_odezero.h"
#include "cr3bp_derivatives_rel.h"

//GSL
#include "../lib/gsl/gsl_odeiv2.h"
#include "../lib/gsl/gsl_errno.h"

//-------------------------------------------------------------------------
// C function
//-------------------------------------------------------------------------
void ode78_bcp_event_rel(double *t,                    //current time
                         double *y,                    //current state
                         double **ye,                  //event states
                         double *te,                   //event times
                         double const *y0,             //initial condition
                         double const t0,              //initial time
                         double tf,                    //final time
                         int nvar,                     //number of state variables
                         struct value_params *val_par, //event structure
                         int *nEvents,                 //the effective number of events found
                         double *param);               //bcp parameters


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
// 8. struct struct_event the event structure 
//-------------------------------------------------------------------------
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]);
