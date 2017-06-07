% ode78_cr3bp.m Help file for ode78_cr3bp MEX-file.
%  ode78_cr3bp.c - a Runge-Kutta 7/8 integrator of the CR3BP vector field.
% 
%  [TE, YE] = ODE78_CR3BP(TSPAN, Y0, MU) integrates the equations of
%  motion of the the CR3BP with a mass ratio MU, from the initial 
%  conditions  Y0, on the interval TSPAN = [T0 TF]. If Y0 is of size 6,
%  only the state is integrated. If Y0 is of size 6+36, both the state and
%  the state transition matrix are integrated. The final time is stored in
%  TE = TF, and the final state in YE.
%
%  [TE, YE, T, Y] = ODE78_CR3BP(TSPAN, Y0, MU) does the same, but the state
%  and time are stored along the trajectory on a variable grid and provided
%  to the user via the outputs T and Y.
% 
%  The input and output states are given in CR3BP synodical 
%  reference frame as well as CR3BP normalized units. For more
%  details, see Koon et al. 2006, chapter 2, <a href="matlab: 
%  web('http://www.cds.caltech.edu/~marsden/volume/missiondesign/KoLoMaRo_DMissionBook_2011-04-25.pdf','-browser')">(link)</a>.
%
%  BLB 2015