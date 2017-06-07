% ode78_cr3bp_event.m Help file for ode78_cr3bp_event MEX-file.
%  ode78_cr3bp_event.c - a Runge-Kutta 7/8 integrator of the CR3BP 
%  vector field, with event handling.
% 
%  [TE, YE] = ODE78_CR3BP_EVENT(TSPAN, Y0, MU, EVENT) integrates the 
%  equations of motion of the the CR3BP with a mass ratio MU, from the
%  initial conditions  Y0, on the interval TSPAN = [T0 TF]. If Y0 is of 
%  size 6, only the state is integrated. If Y0 is of size 6+36, both the 
%  state and the state transition matrix are integrated. During the 
%  integration, the routine looks for specific event occurences, 
%  as defined in the structure EVENT, which should be of the form given
%  by the routine INIT_EVENT.
%  Since an example is worth a thousand words, here is what the EVENT
%  structure looks like when the user wants to detect all the intersections
%  with the plane x = xc.
%
%   - EVENT.type  = 'X_SECTION' = cst.event.type.X_SECTION
%   - EVENT.value = xc
%   - EVENT.dim   = 1 (X dimension, redundancy with EVENT.type)
%   - EVENT.direction   = -1 (decreasing crossings: x = xc with dot(x) < 0)
%                          0 (all crossings)
%                         +1 (increasing crossings: x = xc with dot(x) > 0) 
%
%   See and use INIT_EVENT routine to properly initialize theses
%   structures.
%   Note that the field EVENT.isterminal is NOT used by ODE78_CR3BP_EVENT.
%
%   At the end of the integration, the times and states corresponding to
%   the occurences of the event are stored in [TE, YE].
%
%  [TE, YE, T, Y] = ODE78_CR3BP_EVENT(TSPAN, Y0, MU, EVENT) does the same,
%  but the state and time are stored along the trajectory on a variable 
%  grid and provided to the user via the outputs T and Y.
%
%  WARNING: the order of the outputs is reversed with respect to the usual
%  outputs from MATLAB ode routines used with event.
% 
%  The input and output states are given in CR3BP synodical 
%  reference frame as well as CR3BP normalized units. For more
%  details, see Koon et al. 2006, chapter 2, <a href="matlab: 
%  web('http://www.cds.caltech.edu/~marsden/volume/missiondesign/KoLoMaRo_DMissionBook_2011-04-25.pdf','-browser')">(link)</a>.
%
%  See also INIT_EVENT, ODE78_CR3BP, MANIFOLD_BRANCH_COMPUTATION, EVENT
%
%  BLB 2015