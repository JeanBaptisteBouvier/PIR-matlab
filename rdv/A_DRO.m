function [dro, default, cr3bp] = A_DRO(Ax, default, cst)
% Ax between 16'000 km and 349'000km of radius

%% Inner changes from default parameters

% 1. Choose between MEX and MATLAB
%--------------------------------------------------------------------------
default.computation.type = cst.computation.MEX;

% 3. Decomment the next lines to change the plotting options
%--------------------------------------------------------------------------
default.plot.XZ             = false; % plot also the results in X-Z plane
default.plot.YZ             = false; % plot also the results in Y-Z plane
default.plot.diff_corr      = false; % plot the differential correction steps
default.plot.orbit          = false; % prevent the plotting
default.plot.firstPrimDisp  = true;  % display the Earth

% 4. See parameters_default_init.m to see other options
%--------------------------------------------------------------------------

%% Environment init
cr3bp = init_CR3BP('EARTH', 'MOON', default);
load Dro_data;

%% Orbit init & computation for a planar lyapunov orbit

x = 1 - cr3bp.mu - (Ax/cr3bp.L);


if x < 0.077868572624673 || x > 0.948773910611288 % extrem values of Dro_data.x
    error('Ax value is out of range');         
end

dro = init_orbit(cr3bp, ...        % Parent CR3BP
    cr3bp.l2, ...                  % Parent libration point
    cst.orbit.type.DRO, ...        % Based on a Planar lyapunov orbit
    cst.orbit.family.PLANAR, ...   % Planar class (useless here, since it is a planar lyapunov orbit
    Ax, ...                        % Of Ax extension ~ 80 000 km
    cst);                          % Numerical constants

C = interp1(Dro_data.x, Dro_data.C, x);
v_y = sqrt(abs( x^2  - C + 2*( (1-cr3bp.mu)/abs(x+cr3bp.mu) + cr3bp.mu/abs(x-1+cr3bp.mu) ) ));
yvg = [ x; 0; 0; 0 ; v_y ; 0];
dro = orbit_refinement(cr3bp, dro, default, yvg, cst);




%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
    figure(4);
    view([-47 28]);
end

end
