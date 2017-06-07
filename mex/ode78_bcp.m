%  ode78_bcp.m Help file for ode78_bcp MEX-file.
%  ode78_bcp.c - a Runge-Kutta 7/8 integrator of the Bicircular
%  Sun-Earth-Moon problem.
% 
%  [TE, YE] = ODE78_BCP(TSPAN, Y0, MU, THETA0, MS, AS, OMEGAS)
%  integrates the equations of motion of the the BICIRCULAR SUN-EARTH-MOON 
%  problem with:
%           - MU the Earth-Moon mass ratio,
%           - MS the mass of the Sun,               all in Earth-Moon units
%           - AS the semi-major axis of the Sun,
%     	    - OMEGAS the mean motion of the Sun,
%  from the initial conditions  Y0, on the interval TSPAN = [T0 TF]. 
%  If Y0 is of size 6, only the state is integrated. If Y0 is of size 6+36, 
%  both the state and the state transition matrix are integrated. 
%  The final time is stored in TE = TF, and the final state in YE.
%
%  [TE, YE, T, Y] = ODE78_BCP(TSPAN, Y0, MU, THETA0, MS, AS, OMEGAS) does 
%  the same, but the state and time are stored along the trajectory on a 
%  variable grid and provided to the user via the outputs T and Y.
%
%  The input and output states are given in Earth-Moon synodical 
%  reference frame as well as Earth-Moon normalized units. For more
%  details, see Koon et al. 2006, chapter 2, <a href="matlab: 
%  web('http://www.cds.caltech.edu/~marsden/volume/missiondesign/KoLoMaRo_DMissionBook_2011-04-25.pdf','-browser')">(link)</a>.
% 
%  BLB 2015