/*=========================================================================
 * ode78_bcp_event.c mex file to generate
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
#include "ode78.h"
#include "custom_ode.h"
#include "custom_odezero.h"
#include "cr3bp_derivatives.h"

//GSL
#include "../lib/gsl/gsl_odeiv2.h"
#include "../lib/gsl/gsl_errno.h"

//-------------------------------------------------------------------------
// C function
//-------------------------------------------------------------------------
void ode78_bcp_event(double *t,                    //current time
                     double *y,                    //current state
                     double **ye,                  //event states
                     double *te,                   //event times
                     double const *y0,             //initial condition
                     double const t0,              //initial time
                     double tf,                    //final time
                     int nvar,                     //number of state variables
                     struct value_params *val_par, //event structure
                     int *nEvents,                 //the effective number of events found
                     double *param)                //bcp parameters
{
    //---------------------------------------------------------------------
    // Initialize the integration structures
    //---------------------------------------------------------------------
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    
    
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent; //Brent-Dekker root finding method

    //General ode structure
    custom_ode_structure ode_s;
    switch(nvar)
    {
        case 6:
        {
            init_ode_structure(&ode_s,   T, T_root, eps_AbsTol, eps_RelTol, eps_Root, eps_Diff, 6, 1e-6, bcp_derivatives_6, NULL, param);
            break;
        }
        case 42:
        {
            init_ode_structure(&ode_s,   T, T_root, eps_AbsTol, eps_RelTol, eps_Root, eps_Diff, 42, 1e-6, bcp_derivatives_42, NULL, param);
            break;
        }
        default:   //if nvar !=6 && nvar != 42, return without integration
        {
            mexPrintf("ode78_bcp. Error: wrong number nvar of state variables. nvar must 6 or 42. return.");
            return;
        }
    
    }
    
    
    //---------------------------------------------------------------------
    // Initalization of the state
    //---------------------------------------------------------------------
    *t = t0;
    for(int i = 0; i < nvar; i++) y[i] = y0[i];   
    
    //---------------------------------------------------------------------
    // Integration until t = t1, and saving the events
    //---------------------------------------------------------------------
    struct value_function fvalue;
    fvalue.val_par = val_par;
    
    switch(val_par->type)
    {
        case 'X':
        case 'Y':
        case 'Z':
        {
            fvalue.value   = &linear_intersection; 
            break;
        }
        
        case 'A':
        {
            fvalue.value   = &angle_intersection; 
            break;
        }
        
        case 'F':
        {
            fvalue.value   = &null_flight_path_angle; 
            break;
        }
        
        default:
        {
            mexPrintf("Bad event type. return.");
            return;
        }
    }
    
    
    //Starting in the right direction
    ode_s.d->h = (tf>t0) ? fabs(ode_s.d->h) : -fabs(ode_s.d->h);
    ode_s.h = (tf>t0) ? fabs(ode_s.h) : -fabs(ode_s.h);
    
    //Apply custom_odezero_2
    *nEvents = custom_odezero_2(y, ye, te, t0, tf, &ode_s, fvalue);
    
    // custom_odezero_2 return the last indix rather than the number of events
    // So we need to add one event
    *nEvents = *nEvents+1;
}


//-------------------------------------------------------------------------
// The gateway function.
// The input must be, in that order:
// 1. [t0, tf] the time span
// 2. double y0[6 or 42], the initial state
// 3. double mu, the cr3bp mass ratio
// 4. double theta0, the initial phase of the Sun
// 5. double ms, the mass of the Sun
// 6. double as, the semi-major axis of the Sun
// 7. double omS, the mean motion of the Sun
// 8. struct struct_event the event structure 
//-------------------------------------------------------------------------
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //---------------------------------------------------------------------
    // Check for proper number of arguments
    //---------------------------------------------------------------------
    if(nrhs!=8) {
        mexErrMsgIdAndTxt("custom:ode78_bcp_event:nrhs","8 inputs required.");
    }
    if(!(nlhs == 2 || nlhs == 4)) {
        mexErrMsgIdAndTxt("custom:ode78_bcp_event:nlhs",
       "2 or 4 outputs required: 2 if only the final state is desired, 4 if the state on a given time grid is also desired.");
    }
 
    //---------------------------------------------------------------------
    // Retrieve the variables:
    // 1. t0, the initial time
    // 2. tf, the final time
    // 3. y0, the initial state
    // 4. nvar, the number of state variables
    // 5. mu, the cr3bp mass ratio
    // 6. om0, the initial phase angle of the Sun
    // 7. struct struct_event the event structure 
    //---------------------------------------------------------------------
    double *ts = mxGetPr(prhs[0]);
    int nts    = mxGetN(prhs[0]);
    double *y0 = mxGetPr(prhs[1]);
    int nvar   = mxGetN(prhs[1]);
    int mvar   = mxGetM(prhs[1]);
    double mu  = mxGetScalar(prhs[2]);
    double th0 = mxGetScalar(prhs[3]);
    double ms  = mxGetScalar(prhs[4]); 
    double as  = mxGetScalar(prhs[5]);
    double omS = mxGetScalar(prhs[6]);
    
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
        mexPrintf("nvar = %d\n", nvar);
        mexErrMsgIdAndTxt("custom:ode78_bcp:nts","The state vector (second input) must be of size either 6 or 42.");  
    }
    
    //---------------------------------------------------------------------
    // Build the times
    //---------------------------------------------------------------------
    double t0 = ts[0];
    double tf = ts[1];
    
    //---------------------------------------------------------------------
    // particular case of the event structure:
    //---------------------------------------------------------------------
    // Get input arguments
    mwSize nfields = mxGetNumberOfFields(prhs[7]);           //number of fields
    mwSize NStructElems = mxGetNumberOfElements(prhs[7]);    //number of elements (must be 1)
    
    if(NStructElems != 1)
    {
        mexErrMsgIdAndTxt("custom:ode78_bcp_event:NStructElems",
                "Only one element is required in the event structure.");
    }
    
    
    if(nfields != 7)
    {
        printf(" nfields = %d", nfields);
        mexErrMsgIdAndTxt("custom:ode78_bcp_event:nfields",
                "The event structure must contain 7 fields: type, isterminal, max_events, direction, dim, value, and center.");
    }
    
    //---------------------------------------------------------------------
    //Get the different fields in the event structure val_par = prhs[7]:
    // val_par.max_events
    // val_par.direction
    // val_par.dim
    // val_par.value
    //---------------------------------------------------------------------
    mxArray *tmp;
    int max_events, direction, dim;
    double value;
    double *center = (double*) calloc(3, sizeof(double));
    

    //0. Type
    tmp = mxGetField(prhs[7], 0, "type");
    int tmpN = mxGetN(tmp)+1;
    char type[tmpN];
    mxGetString(tmp, type, tmpN);
    //mexPrintf("type = %s\n", type);
    
    //Check that the event type is not free.
    if (strcmp(type, "FREE") == 0) {
        mexErrMsgIdAndTxt("custom:ode78_bcp_event:type",
                "The event structure as type FREE. Use ode78_bcp instead.");
    }
    
    //1. val_par.max_events
    tmp = mxGetField(prhs[7], 0, "max_events");
    max_events = (int) mxGetScalar(tmp);
    //mexPrintf("max_events = %d\n", max_events); 
    
    //2. val_par.direction
    tmp = mxGetField(prhs[7], 0, "direction");
    direction = (int) mxGetScalar(tmp);
    //mexPrintf("direction = %d\n", direction); 
    
    //3. val_par.dim
    tmp = mxGetField(prhs[7], 0, "dim");
    dim = (int) mxGetScalar(tmp)-1;
    //mexPrintf("dim = %d\n", dim); 
    
    //4. val_par.value
    tmp = mxGetField(prhs[7], 0, "value");
    value = mxGetScalar(tmp);
    //mexPrintf("value = %5.5f\n", value); 
    
    //5. val_par.center
    tmp = mxGetField(prhs[7], 0, "center");
    center = mxGetPr(tmp);
    //mexPrintf("center = (%5.5f, %5.5f, %5.5f)\n", center[0], center[1], center[2]);
    

    //---------------------------------------------------------------------
    // Create the corresponding C structure
    //---------------------------------------------------------------------
    struct value_params val_par;
    val_par.max_events = max_events;
    val_par.direction  = direction;
    val_par.dim        = dim;
    val_par.value      = value;
    val_par.center     = center;
    //Note: only the first letter is used! 
    //So do not create an event with the same initials as another one...
    val_par.type       = (int) *type; 
    
    //---------------------------------------------------------------------
    // Number of events
    //---------------------------------------------------------------------
    int MAX_EVENTS = 1000;  //maximum allowed is 1000
    int nEv; //actual number of events stored at the end of the computation
    
    if(max_events > MAX_EVENTS)
    {
        mexErrMsgIdAndTxt("custom:ode78_cr3bp_events:maxEv","Number of required events must be < %d.", MAX_EVENTS);
    }
    
    
    //---------------------------------------------------------------------
    // Initialization
    //---------------------------------------------------------------------
    //Current state & time
    double t, y[nvar];
    //Event states
    double **ye = (double **) calloc(nvar, sizeof(double*));
    for(int n = 0; n < nvar; n++) ye[n] = (double*) calloc(MAX_EVENTS+1, sizeof(double));
    //Event time
    double *te = (double*) calloc(MAX_EVENTS+1, sizeof(double));
    
    //---------------------------------------------------------------------
    // BCP params
    //---------------------------------------------------------------------
    //Sun parameters
    double *params = (double*) calloc(5, sizeof(double));
    params[0] = mu;
    params[1] = th0;
    params[2] = ms;
    params[3] = as;
    params[4] = omS;
    
    //---------------------------------------------------------------------
    // Integration
    //---------------------------------------------------------------------
    ode78_bcp_event(&t, y, ye, te, y0, t0, tf, nvar, &val_par, &nEv, params);
     
    //print the number of events found
    //mexPrintf("ode78_bcp_event. %d were found.\n", nEv);
    
    //---------------------------------------------------------------------
    // If it is desired, store the state on a time grid up to te[end]
    //---------------------------------------------------------------------
    int nGrid = 1000, nV;
    double **yv, *tv;
    //New final time
    tf = te[nEv-1];
    //Integration, if necessary
    if(nlhs == 4)
    {
        //Old version with fixed grid
        //ode78_cr3bp_vec(&t, y, tv, yv, nGrid, y0, t0, tf, nvar, &mu);
        
        //New version
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
            nV = ode78_bcp_vec_var(&t, y, tv, yv, nGrid, y0, t0, tf, nvar, params);
            
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
    
    //---------------------------------------------------------------------
    // Output
    //---------------------------------------------------------------------
    // create the output matrices
    plhs[0] = mxCreateDoubleMatrix(nEv, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nEv, nvar, mxREAL);
    
    //get a pointer to the real data in the output matrix
    double *tout = mxGetPr(plhs[0]);
    double *yout = mxGetPr(plhs[1]);
      
    //Store the event states
    int indix = 0;
    for(int i = 0; i < nvar; i++)
    {
        for(int k = 0; k < nEv; k++) 
        {
            yout[indix++] = ye[i][k];  
        }       
    }
    for(int k = 0; k < nEv; k++) tout[k] = te[k]; 
    
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

        //Print the gap between the last event and the grid computation
        //         mexPrintf("ode78_bcp_event. Discrepancy between last event and the grid computation\n");
        //         mexPrintf("ye    yvec    gap\n");
        //         for(int i = 0; i < nvar; i++)
        //         {
        //             mexPrintf("%5.5f     %5.5f      %5.5e \n", ye[i][nEv-1], yv[i][nGrid], ye[i][nEv-1]-yv[i][nGrid]);
        //         }
    }
    
    //---------------------------------------------------------------------
    // Free memory
    //---------------------------------------------------------------------
    if(nlhs == 4)
    {
        free(tv);
        for (int i=0; i<=nvar; i++) free(yv[i]); free(yv);
    }
    free(te);
    free(params);
    for (int i=0; i<=nvar; i++) free(ye[i]); free(ye);
}