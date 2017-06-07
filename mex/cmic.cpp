/*=========================================================================
 * cm.cpp 
 *
 * Read a given manifold in semi-analytical form
 * 
 * author:  BLB
 * year:    2016
 * version: 1.0
 *=======================================================================*/

//-------------------------------------------------------------------------
// Headers
//-------------------------------------------------------------------------
//Mex
#include "mex.h"


//Custom
#include "Cpp/FTA.h"
#include "Cpp/Ofsc.h"
#include "Cpp/Oftsc.h"
#include "Cpp/env.h"
#include "Cpp/pmcoc.h"

//-------------------------------------------------------------------------
// Available coordinates system
//-------------------------------------------------------------------------
#define NC 0
#define SYS 1
#define VSYS 2

//-------------------------------------------------------------------------
// Compute the intial conditions (IC) z1NC inside the center manifold
// \param st0   : the IC in RCM coordinates
// \param  t0   : the initial time
// \param z1NC : the IC output in NC coordinates
//-------------------------------------------------------------------------
void computeCM(double *st0, double t0, double  *z1NC)
{ 
    //------------------------------------------
    // Initialize the manifold
    //------------------------------------------
    vector<Oftsc> CM(6);
    readVOFTS_bin(CM, SEML.cs.F_PMS+"W/W", OFS_ORDER);
    
    //------------------------------------------
    // RCM to NC for NC initial conditions
    //------------------------------------------
    Ofsc AUX;            //OFS temp variable
    // RCM to NC
    RCMtoNC(st0, t0, SEML.us.n, OFTS_ORDER, OFS_ORDER, CM, AUX, z1NC, false);
}

//-------------------------------------------------------------------------
// The gateway function.
//
// Compute the initial conditions (IC) z1sys inside the center manifold
// in a given coordinates system. Either:
//  - NC: Normalized-Centered coordinates (x, p) (canonical)
//  - SYS: Three-Body Problem (EM or SEM) coordinates (x, p) (canonical) with the US convention (the Earth has a negative abscissa)
//  - VSYS: Three-Body Problem (EM or SEM) coordinates (x, v) (non canonical) again with the US convention (the Earth has a negative abscissa)
//
//-------------------------------------------------------------------------
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //---------------------------------------------------------------------
    // Set the global variables
    //---------------------------------------------------------------------
    MODEL_TYPE = 0;
    OFS_ORDER  = 0;
    REDUCED_NV = 4;
    
    //---------------------------------------------------------------------
    // Check for proper number of arguments
    //---------------------------------------------------------------------
    if(nrhs!=5) {
        mexErrMsgIdAndTxt("custom:cm:nrhs","5 input required.");
    }
    
    
    //---------------------------------------------------------------------
    // Retrieve the variables:
    // 1. s0 the initial configuration
    // 2. nv number of elements in s0
    //---------------------------------------------------------------------
    double *st0      = mxGetPr(prhs[0]);
    int nvar         = mxGetN(prhs[0]);
    double t0        = mxGetScalar(prhs[1]);
    int order        = (int) mxGetScalar(prhs[2]);
    int li_EM        = (int) mxGetScalar(prhs[3]);
    int outputType   = (int) mxGetScalar(prhs[4]);
    
    //---------------------------------------------------------------------
    // Checks
    //---------------------------------------------------------------------
    //Check that the number of elements in s0 is 4 
    //(dimension of the center manifold)
    if(nvar!=REDUCED_NV) {
        mexErrMsgIdAndTxt("custom:cm:nvar","nvar!=REDUCED_NV.");
    }
    
    //Check that the order if smaller than the maximum order
    //allowed
    if(order>40) {
        mexErrMsgIdAndTxt("custom:cm:order","order>20.");
    }
    
    //Check that the outputType is valid
    if(outputType > 3){
        mexErrMsgIdAndTxt("custom:cm:outputType","outputType>3.");
    }
    
    //---------------------------------------------------------------------
    // Set OFTS_ORDER from order
    //---------------------------------------------------------------------
    OFTS_ORDER = order;
        
    //---------------------------------------------------------------------
    // Retrieve the configuration parameters
    //---------------------------------------------------------------------
    int compType  = 3;
    int model     = 0;
    int dcs       = 0;
    int isNorm    = 1;
    int li_SEM    = 2;
    int pms       = 0;
    int mType_EM  = 0;
    int mType_SEM = 0;
    int storage   = 0;

    //---------------------------------------------------------------------
    // Computation
    //---------------------------------------------------------------------

    //------------------------------------------
    // Initialization of the environnement
    // Mandatory to perform any computation except qbtbp(int)
    //------------------------------------------
    init_env(li_EM, li_SEM, isNorm, model, dcs, pms, mType_EM, mType_SEM);


    //------------------------------------------
    // Computation of the initial conditions
    //------------------------------------------
    double  z1NC[6];    //NC
    computeCM(st0, t0, z1NC);
    
    //------------------------------------------
    // Postprocess depending on the value outputType
    //------------------------------------------
    double  z1sys[6];
    
    switch(outputType)
    {
        case NC:
            for(int i = 0; i < 6; i++) z1sys[i] = z1NC[i];
            break;
        case SYS:
            // NC to SYS
            NCtoSYS(t0, z1NC, z1sys, &SEML);
            
            // Inverse some elements to get the good signs
            z1sys[0] = -z1sys[0];
            z1sys[1] = -z1sys[1];
            z1sys[3] = -z1sys[3];
            z1sys[4] = -z1sys[4];
    
    
            break;
        case VSYS:
            // NC to SYS
            NCtoSYS(t0, z1NC, z1sys, &SEML);
            
            // Inverse some elements to get the good signs
            z1sys[0] = -z1sys[0];
            z1sys[1] = -z1sys[1];
            z1sys[3] = -z1sys[3];
            z1sys[4] = -z1sys[4];
            
            //From momenta to velocities
            z1sys[3] = z1sys[3] + z1sys[1]; //vx = px + y
            z1sys[4] = z1sys[4] - z1sys[0]; //vy = py - x
    }
        
    
    //---------------------------------------------------------------------
    // Output
    //---------------------------------------------------------------------
    //Create the output matrices
    plhs[0] = mxCreateDoubleMatrix(1, 6, mxREAL);
    //Get a pointer to the real data in the output
    double *yout = mxGetPr(plhs[0]);
    //Store the initial state
    for(int i = 0; i < 6; i++) yout[i] = z1sys[i];
}