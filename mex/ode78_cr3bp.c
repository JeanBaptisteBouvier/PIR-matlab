/*=========================================================================
 * ode78_cr3bp.c mex file to generate
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

//Custom
#include "C/ode78.h"


//-------------------------------------------------------------------------
// The gateway function.
// The input must bet, in that order:
// 1. [t0, tf] the time span
// 2. double y0[6 or 42], the initial state
// 3. double mu, the cr3bp mass ratio
//-------------------------------------------------------------------------
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    //---------------------------------------------------------------------
    // Check for proper number of arguments
    //---------------------------------------------------------------------
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("custom:ode78_cr3bp:nrhs","3 inputs required.");
    }
    if(!(nlhs == 2 || nlhs == 4)) {
        mexErrMsgIdAndTxt("custom:ode78_cr3bp:nlhs",
                "2 or 4 outputs required: 2 if only the final state is desired, 4 if the state on a given time grid is also desired.");
    }
    
    //---------------------------------------------------------------------
    // Retrieve the variables:
    // 1. [t0, tf] the time span
    // 2. double y0[6 or 42], the initial state
    // 3. double mu, the cr3bp mass ratio
    //---------------------------------------------------------------------
    double *ts = mxGetPr(prhs[0]);
    int nts    = mxGetN(prhs[0]);
    double *y0 = mxGetPr(prhs[1]);
    int nvar   = mxGetN(prhs[1]);
    int mvar   = mxGetM(prhs[1]);
    double mu  = mxGetScalar(prhs[2]);
    
    //---------------------------------------------------------------------
    // Set the size of the state as the max(nvar, mvar);
    //---------------------------------------------------------------------
    nvar = nvar > mvar? nvar:mvar;
    
    //---------------------------------------------------------------------
    // Do some checks on the inputs
    //---------------------------------------------------------------------
    if(nts!=2) {
        mexErrMsgIdAndTxt("custom:ode78_bcp:nts","The time vector (first input) must be of size 2: [t0 tf]");
    }
    
    if(nvar!=6 && nvar!=42)
    {
        mexErrMsgIdAndTxt("custom:ode78_bcp:nts","The state vector (second input) must be of size either 6 or 42.");
    }
    
    //---------------------------------------------------------------------
    // Build the times
    //---------------------------------------------------------------------
    double t0 = ts[0];
    double tf = ts[1];
    
    //---------------------------------------------------------------------
    // Integration, for 2 outputs
    //---------------------------------------------------------------------
    double **yv, *tv, t, y[nvar];
    int nGrid = 1000, nV;
    if(nlhs == 2)
    {
        ode78_cr3bp(&t, y, y0, t0, tf, nvar, &mu);
    }
    else if(nlhs == 4)
    {
        //-----------------------------------------------------------------
        // Integration, for 4 outputs
        //-----------------------------------------------------------------        
        do{   
            //-------------------------------------------------------------
            // State will be stored on a given grid, if necessary
            //-------------------------------------------------------------
            //Event states
            yv = (double **) calloc(nvar, sizeof(double*));
            for(int n = 0; n < nvar; n++) yv[n] = (double*) calloc(nGrid+1, sizeof(double));
            //Event time
            tv = (double*) calloc(nGrid+1, sizeof(double));
            
            //-------------------------------------------------------------
            // Integration
            //-------------------------------------------------------------
            nV = ode78_cr3bp_vec_var(&t, y, tv, yv, nGrid, y0, t0, tf, nvar, &mu);
            
            //-------------------------------------------------------------
            // Check that the grid is not overflowed
            //-------------------------------------------------------------
            if(nV < nGrid)
            {
                nGrid = nV;   
                break;
            }
            else
            {
                free(tv);
                for (int i=0; i<=nvar; i++) free(yv[i]); free(yv);
                nGrid *= 2;
            }
            
        }while(nGrid < 5e6);
    }
    
    if(nGrid >= 5e6)
    {
        mexErrMsgIdAndTxt("custom:ode78_cr3bp_events:maxCapacity","Maximum capacity reached. Try with a smaller integration time.");
    }
    
    //Old version with fixed grid
    //if(nlhs == 4) ode78_cr3bp_vec(&t, y, tv, yv, nGrid, y0, t0, tf, nvar, &mu);
    
    //---------------------------------------------------------------------
    // Output: the final state
    //---------------------------------------------------------------------
    //Create the output matrices
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);           //t
    plhs[1] = mxCreateDoubleMatrix(1, nvar, mxREAL);        //y
    
    
    //Get a pointer to the real data in the output
    double *tout  = mxGetPr(plhs[0]);
    double *yout  = mxGetPr(plhs[1]);
    
    //Store the final state
    for(int i = 0; i < nvar; i++) yout[i] = y[i];
    tout[0] = t;
    
    //---------------------------------------------------------------------
    // Output: the state on a time grid, if necessary
    //---------------------------------------------------------------------
    if(nlhs == 4)
    {
        //Create the output matrices
        plhs[2] = mxCreateDoubleMatrix(nGrid+1, 1, mxREAL);     //tv
        plhs[3] = mxCreateDoubleMatrix(nGrid+1, nvar, mxREAL);  //yv
        
        //Get a pointer to the real data in the output
        double *tvout = mxGetPr(plhs[2]);
        double *yvout = mxGetPr(plhs[3]);
        
        //Store the state on the grid [0,..., nGrid]
        int indix = 0;
        for(int i = 0; i < nvar; i++)
        {
            for(int k = 0; k <= nGrid; k++)
            {
                yvout[indix++] = yv[i][k];
            }
        }
        for(int k = 0; k <= nGrid; k++) tvout[k] = tv[k];
        
    }
    
    //---------------------------------------------------------------------
    // Free memory
    //---------------------------------------------------------------------
    if(nlhs == 4)
    {
        free(tv);
        for (int i=0; i<=nvar; i++) free(yv[i]); free(yv);
    }
}