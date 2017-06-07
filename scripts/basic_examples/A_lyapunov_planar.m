%--------------------------------------------------------------------------
% This matlab file allows to build:
% - an EML2 planar lyapunov orbit (~8000 km of radius)
%
% BLB 2017
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Inner changes from default parameters

% 1. Decomment the next line to use only MATLAB routines (very slow!)
%--------------------------------------------------------------------------
default.computation.type = cst.computation.MATLAB;

% 3. Decomment the next lines to change the plotting options
%--------------------------------------------------------------------------
default.plot.XZ             = false; % plot also the results in X-Z plane
default.plot.YZ             = false; % plot also the results in Y-Z plane
default.plot.diff_corr      = false; % plot the differential correction steps

% 4. See parameters_default_init.m to see other options
%--------------------------------------------------------------------------

%% Environment init
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%% Orbit init & computation for a planar lyapunov orbit


%Computation choice 1 : another DRO
Ax = 3000 ;

plyap = init_orbit(cr3bp, ...      % Parent CR3BP
    cr3bp.l2, ...                  % Parent libration point
    cst.orbit.type.PLYAP, ...      % Planar lyapunov orbit
    cst.orbit.family.PLANAR, ...   % Planar class (useless here, since it is a planar lyapunov orbit
    Ax, ...                        % Of Ax "radius"
    cst);                          % Numerical constants

plyap = orbit_computation(cr3bp, plyap, default, cst); %calcul of yvg and use orbit_refinement



%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
    figure(4);
    view([-47 28]);
end
