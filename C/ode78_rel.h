/*=========================================================================
 * ode78_rel.h 
 * 
 * C header file that contains the basic routines that appear in the mex files
 *  - ode78_cr3bp_rel.c
 *  - ode78_cr3bp_event_rel.c
 *
 * author:  BLB
 * year:    2015
 * version: 1.0
 * author: JBB
 * year:   2017
 * version: 2.0
 *=======================================================================*/

//-------------------------------------------------------------------------
// Headers
//-------------------------------------------------------------------------
//Mex
#include "mex.h"

//Custom
#include "custom_ode.h"
#include "cr3bp_derivatives_rel.h"

//GSL
#include "../lib/gsl/gsl_odeiv2.h"
#include "../lib/gsl/gsl_errno.h"

//Precisions
extern double eps_AbsTol;
extern double eps_RelTol;
extern double eps_Diff;
extern double eps_Root;

/*-------------------------------------------------------------------------
 * Integration of the CRTBP vector field of mass ratio mu, from the initial
 * state y0(t0) to the final time tf. 
 * 
 * The inputs are:
 *  - double const *y0,  the initial conditions
 *  - double const t0,   initial time
 *  - double tf,         the final time
 *  - int nvar,          the number of state variables
 *  - double *mu,        the cr3bp mass ratio
 *
 * The outputs are:
 *  - double *t,  the current time (should be equal to tf at the end).
 *  - double *y,  the current state, equal to y(tf) at the end of the process
 *
 *  Remarks:
 *   - The number of variables is either 12 (relative state & target state)
 *          or 48 (full state + State Transition Matrix)
 *
 *   - The vector field is computed through either the routine 
 *      cr3bp_derivatives_12 or cr3bp_derivatives_48
 * ----------------------------------------------------------------------*/
void ode78_cr3bp_rel(double *t,                  //current time
                     double *y,                  //current state
                     double const *y0,           //initial condition
                     double const t0,            //initial time
                     double tf,                  //final time
                     int nvar,                   //number of state variables
                     double *mu);                //cr3bp mass ratio


/*-------------------------------------------------------------------------
 * Integration of the CRTBP vector field of mass ratio mu, from the initial
 * state y0(t0) to the final time tf. The state is stored along the
 * trajectory on a given grid.
 * 
 * The inputs are:
 *  - int nGrid,         nGrid+1 is the number of point on the grid
 *  - double const *y0,  the initial conditions
 *  - double const t0,   initial time
 *  - double tf,         the final time
 *  - int nvar,          the number of state variables
 *  - double *mu,        the cr3bp mass ratio
 *
 * The outputs are:
 *  - double *t,  the current time (should be equal to tf at the end).
 *  - double *y,  the current state, equal to y(tf) at the end of the process
 *  - double *tv, the time on the grid [0, ..., nGrid]
 *  - double *yv, the state on the grid [0, ..., nGrid]
 *
 *  Remarks:
 *   - The number of variables is either 6 (full state) 
 *          or 48 (full state + State Transition Matrix)
 *
 *   - The vector field is computed through either the routine 
 *      cr3bp_derivatives_12 or cr3bp_derivatives_48
 * ----------------------------------------------------------------------*/
void ode78_cr3bp_vec_rel(double *t,                  //current time
                         double *y,                  //current state
                         double *tv,                 //time on a given grid of size [nGrid+1]
                         double **yv,                //state on a given grid of size [nGrid+1]
                         int    nGrid,               //size of the grid
                         double const *y0,           //initial condition
                         double const t0,            //initial time
                         double tf,                  //final time
                         int nvar,                   //number of state variables
                         double *mu);                //cr3bp mass ratio

/*-------------------------------------------------------------------------
 * Integration of the CRTBP vector field of mass ratio mu, from the initial
 * state y0(t0) to the final time tf. The state is stored along the
 * trajectory on a variable grid.
 * 
 * The inputs are:
 *  - int nGrid,         nGrid+1 is the max number of point on the grid
 *  - double const *y0,  the initial conditions
 *  - double const t0,   initial time
 *  - double tf,         the final time
 *  - int nvar,          the number of state variables
 *  - double *mu,        the cr3bp mass ratio
 *
 * The outputs are:
 *  - double *t,  the current time (should be equal to tf at the end).
 *  - double *y,  the current state, equal to y(tf) at the end of the process
 *  - double *tv, the time on the grid [0, ..., nGrid]
 *  - double *yv, the state on the grid [0, ..., nGrid]
 *  - the function returns the last updated position on the grid (<= nGrid)
 *
 *  Remarks:
 *   - The number of variables is either 12 (relative state & target state) 
 *          or 48 (full state + State Transition Matrix)
 *
 *   - The vector field is computed through either the routine 
 *      cr3bp_derivatives_12 or cr3bp_derivatives_48
 * ----------------------------------------------------------------------*/
int ode78_cr3bp_vec_var_rel(double *t,                  //current time
                            double *y,                  //current state
                            double *tv,                 //time on a given grid of size [nGrid+1]
                            double **yv,                //state on a given grid of size [nGrid+1]
                            int    nGrid,               //size of the grid
                            double const *y0,           //initial condition
                            double const t0,            //initial time
                            double tf,                  //final time
                            int nvar,                   //number of state variables
                            double *mu);                //cr3bp mass ratio

/*-------------------------------------------------------------------------
 * Integration of the BCP vector field of mass ratio mu, from the initial
 * state y0(t0) to the final time tf. 
 * 
 * The inputs are:
 *  - double const *y0,  the initial conditions
 *  - double const t0,   initial time
 *  - double tf,         the final time
 *  - int nvar,          the number of state variables
 *  - double *mu,        the cr3bp mass ratio
 *
 * The outputs are:
 *  - double *t,  the current time (should be equal to tf at the end).
 *  - double *y,  the current state, equal to y(tf) at the end of the process
 *
 *  Remarks:
 *   - The number of variables is either 12 (relative state & target state)
 *          or 48 (full state + State Transition Matrix)
 *
 *   - The vector field is computed through either the routine 
 *      bcp_derivatives_12 or bcp_derivatives_48
 * ----------------------------------------------------------------------*/
void ode78_bcp_rel(double *t,                    //current time
                   double *y,                  //current state
                   double const *y0,           //initial condition
                   double const t0,            //initial time
                   double tf,                  //final time
                   int nvar,                   //number of state variables
                   double *param);             //bcp parameters

/*-------------------------------------------------------------------------
 * Integration of the BCP vector field of mass ratio mu, from the initial
 * state y0(t0) to the final time tf. The state is stored along the
 * trajectory on a given grid.
 * 
 * The inputs are:
 *  - int nGrid,         nGrid+1 is the number of point on the grid
 *  - double const *y0,  the initial conditions
 *  - double const t0,   initial time
 *  - double tf,         the final time
 *  - int nvar,          the number of state variables
 *  - double *mu,        the cr3bp mass ratio
 *
 * The outputs are:
 *  - double *t,  the current time (should be equal to tf at the end).
 *  - double *y,  the current state, equal to y(tf) at the end of the process
 *  - double *tv, the time on the grid [0, ..., nGrid]
 *  - double *yv, the state on the grid [0, ..., nGrid]
 *
 *  Remarks:
 *   - The number of variables is either 12 (relative state & target state) 
 *          or 48 (full state + State Transition Matrix)
 *
 *   - The vector field is computed through either the routine 
 *      bcp_derivatives_12 or bcp_derivatives_48
 * ----------------------------------------------------------------------*/
void ode78_bcp_vec_rel(double *t,                  //current time
                       double *y,                  //current state
                       double *tv,                 //time on a given grid of size [nGrid+1]
                       double **yv,                //state on a given grid of size [nGrid+1]
                       int    nGrid,               //size of the grid
                       double const *y0,           //initial condition
                       double const t0,            //initial time
                       double tf,                  //final time
                       int nvar,                   //number of state variables
                       double *param);             //bcp parameters

/*-------------------------------------------------------------------------
 * Integration of the BCP vector field of mass ratio mu, from the initial
 * state y0(t0) to the final time tf. The state is stored along the
 * trajectory on a variable grid.
 * 
 * The inputs are:
 *  - int nGrid,         nGrid+1 is the max number of point on the grid
 *  - double const *y0,  the initial conditions
 *  - double const t0,   initial time
 *  - double tf,         the final time
 *  - int nvar,          the number of state variables
 *  - double *mu,        the cr3bp mass ratio
 *
 * The outputs are:
 *  - double *t,  the current time (should be equal to tf at the end).
 *  - double *y,  the current state, equal to y(tf) at the end of the process
 *  - double *tv, the time on the grid [0, ..., nGrid]
 *  - double *yv, the state on the grid [0, ..., nGrid]
 *
 *  Remarks:
 *   - The number of variables is either 12 (relative state & target state) 
 *          or 48 (full state + State Transition Matrix)
 *
 *   - The vector field is computed through either the routine 
 *      bcp_derivatives_12 or bcp_derivatives_48
 * ----------------------------------------------------------------------*/
int ode78_bcp_vec_var_rel(double *t,                  //current time
                          double *y,                  //current state
                          double *tv,                 //time on a given grid of size [nGrid+1]
                          double **yv,                //state on a given grid of size [nGrid+1]
                          int    nGrid,               //size of the grid
                          double const *y0,           //initial condition
                          double const t0,            //initial time
                          double tf,                  //final time
                          int nvar,                   //number of state variables
                          double *param);             //bcp parameters