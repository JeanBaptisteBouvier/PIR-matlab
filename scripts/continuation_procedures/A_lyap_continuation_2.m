%--------------------------------------------------------------------------
% Continuation example n°1: Continuation procedure to produce a discrete 
% set within the family of planar lyapunov orbits of EML2
%
% BLB 2015
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Default
default.computation.type = cst.computation.MEX; % MEX files are used to fasten computation

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%Orbit
orbit_1 = init_orbit(cr3bp, cr3bp.l2,  cst.orbit.type.PLYAP, cst.orbit.family.PLANAR, 8400, cst);
orbit_2 = init_orbit(cr3bp, cr3bp.l2,  cst.orbit.type.PLYAP, cst.orbit.family.PLANAR, 8450, cst);

%% Orbit computation
orbit_1 = orbit_computation(cr3bp, orbit_1, default, cst);
orbit_2 = orbit_computation(cr3bp, orbit_2, default, cst);

%% Continuation
maxIter = 10;
%Waitbar
h = waitbar(0,'Computation in progress...');
%Loop
for i = 1:maxIter
    yv = orbit_2.y0 + (orbit_2.y0 - orbit_1.y0);
    orbit_1 = orbit_2;
    orbit_2 = orbit_refinement(cr3bp, orbit_2, default, yv, cst);
    waitbar(i / maxIter);
    period(i) = orbit_2.T;
    abscissa(i) = orbit_2.Axdim;%y0(1);
end
% Close waitbar
close(h)
%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
   figure(4);
   view([-47 28]);
end

%%
figure;
hold on
grid on
plot(abscissa, period, 'o');



%%
default.plot.orbit = false;
options = optimset('Display','iter'); % show iterations
fz = fzero(@(x0)zerof(x0, cr3bp, default, cst), [8400 9000], options);
default.plot.orbit = true;
%%
close all
orbit_2 = init_orbit(cr3bp, cr3bp.l2,  cst.orbit.type.PLYAP, cst.orbit.family.PLANAR, fz, cst);
orbit_2 = orbit_computation(cr3bp, orbit_2, default, cst);

figure;
hold on
grid on
plot(orbit_2.yv(:,1), orbit_2.yv(:,2))
plot(orbit_2.yv(1,1), orbit_2.yv(1,2), 'o')