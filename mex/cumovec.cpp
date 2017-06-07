/*=========================================================================
 * cumovec.cpp
 *
 * Integrates several solutions inside a given center-unstable manifold
 * using the projection method.
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
#include "Cpp/ode.h"
#include "Cpp/vf.h"
#include "Cpp/single_orbit.h"
#include "Cpp/eminsem.h"

//-------------------------------------------------------------------------
// Available coordinates system
//-------------------------------------------------------------------------
#define NC   0
#define SYS  1
#define VSYS 2
#define VSEM 3
//-------------------------------------------------------------------------
//Model
//-------------------------------------------------------------------------
#define M_RTBP  0 // RTBP model indix
#define M_QBCP  1 // QBCP model indix

//Currently not available
//#define M_BCP   2 // BCP model indix
//#define M_ERTBP 3 // ERTBP model indix

//-------------------------------------------------------------------------
//Frameworks
//-------------------------------------------------------------------------
#define F_EM 0
#define F_SEM 1


//-------------------------------------------------------------------------
// The gateway function.
//
// Integrates the initial conditions inside the center manifold using the
// projection method. Inputs:
//  1. st0:   the initial condions in RCM coordinates (4x1 real vector)
//  2. tspan: the desired time span, in the form [t0 t1... tf] (normalized units)
//  3. order: the desired order of the expansions.
//  4. li_EM: the desired libration point (1 or 2)
//  5. outputType: the desired coordinates sytem for the output. Either:
//  - NC: Normalized-Centered coordinates (x, p) (canonical)
//  - SYS: Three-Body Problem (EM or SEM) coordinates (x, p) (canonical) with the US convention (the Earth has a negative abscissa)
//  - VSYS: Three-Body Problem (EM or SEM) coordinates (x, v) (non canonical) again with the US convention (the Earth has a negative abscissa)
//-------------------------------------------------------------------------
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    //---------------------------------------------------------------------
    // Set the global variables
    //---------------------------------------------------------------------
    REDUCED_NV = 5;
    
    //---------------------------------------------------------------------
    // Check for proper number of arguments
    //---------------------------------------------------------------------
    if(nrhs!=12) {
        mexErrMsgIdAndTxt("custom:cm:nrhs","12 input required.");
    }
    
    
    //---------------------------------------------------------------------
    // Retrieve the variables:
    // 1. s0 the initial configuration
    // 2. nv number of elements in s0
    //---------------------------------------------------------------------
    //IC
    double *st0      = mxGetPr(prhs[0]);
    int nvar         = mxGetN(prhs[0]);
    //Time span
    double *tspan    = mxGetPr(prhs[1]);
    int nts          = mxGetN(prhs[1]);
    //Tman
    double tman      = mxGetScalar(prhs[2]);
    //Order
    int order        = (int) mxGetScalar(prhs[3]);
    //Libration point
    int li_EM        = (int) mxGetScalar(prhs[4]);
    //Output Type
    int outputType   = (int) mxGetScalar(prhs[5]);
    //Model
    int model        = (int) mxGetScalar(prhs[6]);
    //Desired framework
    int fwrk         = (int) mxGetScalar(prhs[7]);
    //Man type
    int manType      = (int) mxGetScalar(prhs[8]);
    //kmin
    int kmin         = (int) mxGetScalar(prhs[9]);
    //kmax
    int kmax         = (int) mxGetScalar(prhs[10]);
    //kstep
    int kstep        = (int) mxGetScalar(prhs[11]);
    
    //---------------------------------------------------------------------
    // Checks
    //---------------------------------------------------------------------
    //Check that the number of elements in s0 is 4
    //(dimension of the center manifold)
    if(nvar!=REDUCED_NV) {
        mexErrMsgIdAndTxt("custom:cm:nvar","nvar!=REDUCED_NV.");
    }
    
    //Check that the time span is of size 2
    if(nts!=2) {
        mexErrMsgIdAndTxt("custom:cm:nts","nts!=2.");
    }
    
    //Check that the model is known
    if(model > 1){
        mexErrMsgIdAndTxt("custom:cm:model","model>1.");
    }
    
    int orderMax[2];
    orderMax[0] = 30;   //CRTBP
    orderMax[1] = 20;   //QBCP
    //Check that the order if smaller than the maximum order allowed
    if(order>orderMax[model]) {
        mexErrMsgIdAndTxt("custom:cm:order",
                "orderr is greater than the maximum order currently allowed: 40 for the CRTBP, 20 for the QBCP.");
    }
    
    //Check that the outputType is valid
    if(outputType > 3){
        mexErrMsgIdAndTxt("custom:cm:outputType","outputType>3.");
    }
    
    //Check that the framework is known
    if(fwrk > 1){
        mexErrMsgIdAndTxt("custom:cm:fwrk","fwrk > 1.");
    }
    
    
    
    //---------------------------------------------------------------------
    // Set the global variables
    //---------------------------------------------------------------------
    MODEL_TYPE = model;
    if(MODEL_TYPE == M_RTBP) OFS_ORDER  = 0;
    else OFS_ORDER  = 30;
    
    
    //---------------------------------------------------------------------
    // Set OFTS_ORDER from order
    //---------------------------------------------------------------------
    OFTS_ORDER = order;
    
    
    //---------------------------------------------------------------------
    // Retrieve the configuration parameters
    //---------------------------------------------------------------------
    int isNorm    = 1;
    int li_SEM    = 2;
    int pms       = PMS_GRAPH;
    int mType_EM  = manType;
    int mType_SEM = MAN_CENTER;
    
    //---------------------------------------------------------------------
    // Computation
    //---------------------------------------------------------------------
    
    //---------------------------------------------------------------------
    // Initialization of the environnement
    // Mandatory to perform any computation except qbtbp(int)
    //---------------------------------------------------------------------
    init_env(li_EM, li_SEM, isNorm, model, fwrk, pms, mType_EM, mType_SEM);
    
    //--------------------------------------
    // Center manifolds
    //--------------------------------------
    vector<Oftsc>  CM(6);     ///center manifold in NC coordinates
    vector<Oftsc> CMh(6);     ///center manifold in TFC coordinates
    
    //Read from file
    readVOFTS_bin(CM,  SEML.cs.F_PMS+"W/W",  OFS_ORDER);
    readVOFTS_bin(CMh, SEML.cs.F_PMS+"W/Wh", OFS_ORDER);
    
    //--------------------------------------
    // COC objects
    //--------------------------------------
    matrix<Ofsc>  Mcoc(6,6);    ///COC matrix
    matrix<Ofsc>  Pcoc(6,6);    ///COC matrix (Mcoc = Pcoc*Complex matrix)
    matrix<Ofsc>  MIcoc(6,6);   ///COC matrix = inv(Mcoc)
    matrix<Ofsc>  PIcoc(6,6);   ///COC matrix = inv(Pcoc)
    vector<Ofsc>  Vcoc(6);      ///COC vector
    
    //Read from files
    initCOC(Pcoc, Mcoc, PIcoc, MIcoc, Vcoc, SEML);
    
    //---------------------------------------------------------------------
    // Computation
    //---------------------------------------------------------------------
    //------------------
    // ODE
    //------------------
    OdeStruct driver;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&driver,
            T,               //stepper: int
            T_root,          //stepper: root
            1e-14,           //precision: int_abs
            1e-14,           //precision: int_rel
            1e-13,           //precision: root
            1e-12,           //precision: diffcorr
            6,               //dimension
            1e-6,            //initial int step
            qbfbp_vfn_novar, //vector field
            NULL,            //jacobian
            &SEML);          //current four-body system
    
    //------------------
    //Orbit structure
    //------------------
    Ofsc orbit_ofs(OFS_ORDER);
    SingleOrbit orbit;
    
    //The default interval of projection is set to Tproj = T/5
    double tproj = SEML.us.T/5.0;
    
    //------------------
    // All times (tspan & tman) are multiplied by SEML.us.T
    //------------------
    tspan[0] *= SEML.us.T;
    tspan[1] *= SEML.us.T;
    tman     *= SEML.us.T;
    
    //Init routine
    init_orbit(orbit, &CM, &CMh, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 1, tspan[0], tspan[1], tproj, &driver, &SEML);
    
    //------------------
    //Orbit IC
    //------------------
    orbit_update_ic(orbit, st0, orbit.t0);
    
    //------------------------------------------
    //Plotting parameters
    //------------------------------------------
    double tplot = 1e-2; //we plot every tplot = 0.05
    int nGrid = floor(fabs(orbit.tf - orbit.t0)/tplot);
    double **yNC = dmatrix(0, 5, 0, nGrid);
    double *tv   = dvector(0, nGrid);
    
    //------------------------------------------
    //Integration
    //------------------------------------------
    int status = trajectory_integration_grid(orbit, orbit.t0, orbit.tf, yNC, tv, nGrid, 1);
    
    if(status != 0)
    {
        mexErrMsgIdAndTxt("cmo:precision", "The interval of projection is below the minimum allowed.\n The integration has been stopped and no outputs are given.\n Try a bigger order or a smaller orbit.");
    }
    
    
    //---------------------------------------------------------------------
    // Computation of the manifold leg
    //---------------------------------------------------------------------
    //------------------------------------------
    // Loop parameters
    //------------------------------------------
    int kn    = floor((kmax - kmin)/kstep);
    
    //------------------------------------------
    //Plotting parameters
    //------------------------------------------
    // We need ~ 200 points per Sun period.
    int Nman = 200*floor(fabs(tman)/SEML.us.T);
    
    //to compute the manifold leg
    double **yManNC = dmatrix(0, 5, 0, Nman);
    double *tMan    = dvector(0, Nman);
    
    //to store data
    double ***yManNCStore = d3tensor(0, 5, 0, Nman, 0, kn);
    double **tManStore    = dmatrix(0, Nman, 0, kn);
    int    *ntStore       = ivector(0, kn);
    
    //Current state
    double yc[6], tc;
    
    //RCM state
    double sti[5];
    
    int indix = 0;
    int nt;
    for(int k = kmin; k <= kmax; k+=kstep)
    {
        //----------------------
        // Projection on the center-stable manifold
        //----------------------
        //The current state is the point on the orbit.
        for(int i = 0; i <6; i++) yc[i] = yNC[i][k];
        //The current time is the time on the orbit
        tc = tv[k];
        NCprojCCMtoCUS(yc, tc, sti, orbit, +1e-6, MIcoc, Vcoc);
        
        //----------------------
        // Update the state
        //----------------------
        orbit_update_ic(orbit, sti, tc);
        
        //------------------------------------------
        //Integration
        //------------------------------------------
        nt = trajectory_integration_variable_grid(orbit, tc, tc+tman, yManNC, tMan, Nman, 0);
        //OR
        //trajectory_integration_grid(orbit, tc, tc+tman, yManNC, tMan, Nman, 0);
        //nt = Nman;
        
        //------------------------------------------
        //Store data in NC units
        //------------------------------------------
        for(int i = 0; i <= nt; i++)
        {
            for(int p = 0; p < 6; p++) yManNCStore[p][i][indix] = yManNC[p][i];
            tManStore[i][indix] = tMan[i];
            ntStore[indix] = nt;
        }
        indix++;
    }
    
    
    //------------------------------------------
    // Postprocess depending on the value outputType
    //------------------------------------------
    double **yv     = dmatrix(0, 5, 0, nGrid);
    double ***yvMan = d3tensor(0, 5, 0, Nman, 0, kn);
    double z1NC[6], z1sys[6];
    switch(outputType)
    {
        case NC:
            
            for(int k = 0; k <= nGrid; k++)
            {
                for(int i = 0; i < 6; i++) yv[i][k] = yNC[i][k];
            }
            break;
            
            for(int i = 0; i < 6; i++)
            {
                for(int k = 0; k <= kn; k++)
                {
                    for(int j = 0; j <= ntStore[k]; j++)
                    {
                        yvMan[i][j][k] = yManNCStore[i][j][k];
                    }
                }
            }
            break;
            
        case SYS:
            
            for(int k = 0; k <= nGrid; k++)
            {
                //Update z1NC
                for(int i = 0; i < 6; i++) z1NC[i] = yNC[i][k];
                
                // NC to SYS
                NCtoSYS(tv[k], z1NC, z1sys, &SEML);
                
                // Inverse some elements to get the good signs
                z1sys[0] = -z1sys[0];
                z1sys[1] = -z1sys[1];
                z1sys[3] = -z1sys[3];
                z1sys[4] = -z1sys[4];
                
                //Update yv
                for(int i = 0; i < 6; i++) yv[i][k] = z1sys[i];
            }
            
            //Same for manifold
            for(int k = 0; k <= kn; k++)
            {
                for(int j = 0; j <= ntStore[k]; j++)
                {
                    
                    //Update z1NC
                    for(int i = 0; i < 6; i++)
                    {
                        z1NC[i] = yManNCStore[i][j][k];
                    }
                    
                    // NC to SYS
                    NCtoSYS(tManStore[j][k], z1NC, z1sys, &SEML);
                    
                    // Inverse some elements to get the good signs
                    z1sys[0] = -z1sys[0];
                    z1sys[1] = -z1sys[1];
                    z1sys[3] = -z1sys[3];
                    z1sys[4] = -z1sys[4];
                    
                    //Update yv
                    for(int i = 0; i < 6; i++)
                    {
                        yvMan[i][j][k] = z1sys[i];
                    }
                }
            }
            
            
            break;
            
        case VSYS:
        {
            for(int k = 0; k <= nGrid; k++)
            {
                //Update z1NC
                for(int i = 0; i < 6; i++) z1NC[i] = yNC[i][k];
                
                // NC to SYS
                NCtoSYS(tv[k], z1NC, z1sys, &SEML);
                
                // Inverse some elements to get the good signs
                z1sys[0] = -z1sys[0];
                z1sys[1] = -z1sys[1];
                z1sys[3] = -z1sys[3];
                z1sys[4] = -z1sys[4];
                
                //From momenta to velocities
                z1sys[3] = z1sys[3] + z1sys[1]; //vx = px + y
                z1sys[4] = z1sys[4] - z1sys[0]; //vy = py - x
                
                //Update yv
                for(int i = 0; i < 6; i++) yv[i][k] = z1sys[i];
            }
            
            //Same for manifold
            for(int k = 0; k <= kn; k++)
            {
                for(int j = 0; j <= ntStore[k]; j++)
                {
                    
                    //Update z1NC
                    for(int i = 0; i < 6; i++) z1NC[i] = yManNCStore[i][j][k];
                    
                    // NC to SYS
                    NCtoSYS(tManStore[j][k], z1NC, z1sys, &SEML);
                    
                    // Inverse some elements to get the good signs
                    z1sys[0] = -z1sys[0];
                    z1sys[1] = -z1sys[1];
                    z1sys[3] = -z1sys[3];
                    z1sys[4] = -z1sys[4];
                    
                    //From momenta to velocities
                    z1sys[3] = z1sys[3] + z1sys[1]; //vx = px + y
                    z1sys[4] = z1sys[4] - z1sys[0]; //vy = py - x
                    
                    //Update yv
                    for(int i = 0; i < 6; i++) yvMan[i][j][k] = z1sys[i];
                }
            }
            
            break;
        }
        case VSEM:
        {
            for(int k = 0; k <= nGrid; k++)
            {
                //Update z1NC
                for(int i = 0; i < 6; i++) z1NC[i] = yNC[i][k];
                
                // NC EM to SEM
                if(SEML.fwrk == F_EM) NCEMmtoSEMm(tv[k], z1NC, z1sys, &SEML); //if current fwrk is F_EM, a transformation is needed
                else NCtoSYS(tv[k], z1NC, z1sys, &SEML); //else, we just go to system coordinates
                
                // Inverse some elements to get the good signs
                z1sys[0] = -z1sys[0];
                z1sys[1] = -z1sys[1];
                z1sys[3] = -z1sys[3];
                z1sys[4] = -z1sys[4];
                
                //From momenta to velocities
                z1sys[3] = z1sys[3] + z1sys[1]; //vx = px + y
                z1sys[4] = z1sys[4] - z1sys[0]; //vy = py - x
                
                //Update yv
                for(int i = 0; i < 6; i++) yv[i][k] = z1sys[i];
            }
            
            //Same for manifold
            for(int k = 0; k <= kn; k++)
            {
                for(int j = 0; j <= ntStore[k]; j++)
                {
                    
                    //Update z1NC
                    for(int i = 0; i < 6; i++) z1NC[i] = yManNCStore[i][j][k];
                    
                    // NC EM to SEM
                    if(SEML.fwrk == F_EM) NCEMmtoSEMm(tManStore[j][k], z1NC, z1sys, &SEML); //if current fwrk is F_EM, a transformation is needed
                    else NCtoSYS(tManStore[j][k], z1NC, z1sys, &SEML); //else, we just go to system coordinates
                    
                    // Inverse some elements to get the good signs
                    z1sys[0] = -z1sys[0];
                    z1sys[1] = -z1sys[1];
                    z1sys[3] = -z1sys[3];
                    z1sys[4] = -z1sys[4];
                    
                    //From momenta to velocities
                    z1sys[3] = z1sys[3] + z1sys[1]; //vx = px + y
                    z1sys[4] = z1sys[4] - z1sys[0]; //vy = py - x
                    
                    //Update yv
                    for(int i = 0; i < 6; i++) yvMan[i][j][k] = z1sys[i];
                }
            }
            
            break;
        }
        
    }
    
    
    //---------------------------------------------------------------------
    // Output
    //---------------------------------------------------------------------
    //Create the output matrices
    plhs[0] = mxCreateDoubleMatrix(nGrid+1, 1,   mxREAL);     //tv
    plhs[1] = mxCreateDoubleMatrix(nGrid+1, 6,   mxREAL);     //yv
    plhs[2] = mxCreateDoubleMatrix(Nman+1, kn+1, mxREAL);     //tvMan
    
    //yvMan
    int dims[3];
    dims[0] = Nman+1;
    dims[1] = 6;
    dims[2] = kn+1;
    plhs[3] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    
    //Get a pointer to the real data in the output
    double *tvout = mxGetPr(plhs[0]);
    double *yvout = mxGetPr(plhs[1]);
    double *tvManOut = mxGetPr(plhs[2]);
    double *yvManOut = mxGetPr(plhs[3]);
    
    //Store the state on the grid [0,..., nGrid]
    indix = 0;
    for(int i = 0; i < 6; i++)
    {
        for(int k = 0; k <= nGrid; k++)
        {
            yvout[indix++] = yv[i][k];
        }
    }
    for(int k = 0; k <= nGrid; k++) tvout[k] = tv[k];
    
    //Same for manifold leg
    indix = 0;
    for(int k = 0; k <= kn; k++)
    {
        for(int i = 0; i < 6; i++)
        {
            for(int j = 0; j <= Nman; j++)
            {
                
                yvManOut[indix++] = yvMan[i][j][k];
            }
        }
    }
    
    //Same for manifold leg time
    indix = 0;
    for(int k = 0; k <= kn; k++)
    {
        for(int i = 0; i <= Nman; i++)
        {
            tvManOut[indix++] = tManStore[i][k];
        }
    }
    
    //------------------------------------------
    //Free
    //------------------------------------------
    free_orbit(&orbit);
    free_dmatrix(yNC, 0, 5, 0, nGrid);
    free_dmatrix(yv, 0, 5, 0, nGrid);
    free_dvector(tv, 0, nGrid);
}