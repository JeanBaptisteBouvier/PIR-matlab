%--------------------------------------------------------------------------
% Continuation example n°1: Continuation procedure to produce a discrete 
% set within the family of planar lyapunov orbits of EML2
%
% BLB 2015
%
% JB modifications : not use orbit_computation, directly use
% orbit_refinement so as to choose myself the initialization vector
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;
%JB
default.computation.type = 0;


%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);
dro_data;

%% Orbit 1
k=2;
x_0 = DRO.x_0(k);
Cj = DRO.C_j(k);
Ax = cr3bp.L * (1-cr3bp.mu - x_0);
orbit_1 = init_orbit(cr3bp, cr3bp.l2, cst.orbit.type.PLYAP, cst.orbit.family.PLANAR, Ax, cst);
v_y = sqrt(abs( x_0^2  - Cj + 2*( (1-cr3bp.mu)/abs(x_0+cr3bp.mu) + cr3bp.mu/abs(x_0-1+cr3bp.mu) ) ));
yvg = [ x_0; 0; 0; ;0 ; v_y ; 0];
orbit_1 = orbit_refinement(cr3bp, orbit_1, default, yvg, cst);

%% Orbit 2
k=1; %index of the last dro encountered
x_0 = DRO.x_0(k);
Cj = DRO.C_j(k);
Ax = cr3bp.L * (1-cr3bp.mu - x_0);
orbit_2 = init_orbit(cr3bp, cr3bp.l2, cst.orbit.type.PLYAP, cst.orbit.family.PLANAR, Ax, cst);
v_y = sqrt(abs( x_0^2  - Cj + 2*( (1-cr3bp.mu)/abs(x_0+cr3bp.mu) + cr3bp.mu/abs(x_0-1+cr3bp.mu) ) ));
yvg = [ x_0; 0; 0; ;0 ; v_y ; 0];
orbit_2 = orbit_refinement(cr3bp, orbit_2, default, yvg, cst);





%% Continuation
maxIter = 20;

% Waitbar
h = waitbar(0,'Computation in progress...');

Dro_data_2 = struct('x', [], 'C', [], 'T', []);

% Loop

i = 1;
while abs(orbit_2.y0(1) - (1-cr3bp.mu)) * cr3bp.L > 1800 % don't go to close from the Moon
    yv = orbit_2.y0 + (orbit_2.y0 - orbit_1.y0);
    orbit_1 = orbit_2;
    orbit_2 = orbit_refinement(cr3bp, orbit_2, default, yv, cst);
    waitbar(i / maxIter);
    Dro_data_2.T(i) = orbit_2.T;
    x = orbit_2.y0(1);
    Dro_data_2.x(i) = x;
    Dro_data_2.C(i) = x^2 + 2*( (1-cr3bp.mu)/abs(x+cr3bp.mu) + cr3bp.mu/abs(x-1+cr3bp.mu) ) - orbit_2.y0(5)^2;
    i = 1 + i;
end

% Close waitbar
close(h)

save Dro_data_2 Dro_data_2

%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
   figure(4);
   view([-47 28]);
end

%%
figure;
hold on
grid on
plot(Dro_data_2.x, Dro_data_2.C);

