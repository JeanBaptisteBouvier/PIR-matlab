% 2 hours of calculations needed !!
% Compute the structure Dro_data.mat gathering x, C and T for every DRO
% Dro_data_2 gathers the same data but for closer orbits than the first
% given. It has been obtained with Reverse continuation


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

load Dro_data_2

Dro_data = struct('x', [], 'C', [], 'T', []);
% The first values come from reverse continuation that give smaller orbits
Dro_data.x = fliplr(Dro_data_2.x);
Dro_data.C = fliplr(Dro_data_2.C);
Dro_data.T = fliplr(Dro_data_2.T);

%% Orbit 1
k=1;
x_0 = DRO.x_0(k);
Cj = DRO.C_j(k);
Ax = cr3bp.L * (1-cr3bp.mu - x_0);
orbit_1 = init_orbit(cr3bp, cr3bp.l2, cst.orbit.type.PLYAP, cst.orbit.family.PLANAR, Ax, cst);
v_y = sqrt(abs( x_0^2  - Cj + 2*( (1-cr3bp.mu)/abs(x_0+cr3bp.mu) + cr3bp.mu/abs(x_0-1+cr3bp.mu) ) ));
yvg = [ x_0; 0; 0; ;0 ; v_y ; 0];
orbit_1 = orbit_refinement(cr3bp, orbit_1, default, yvg, cst);

%% Orbit 2
k=2; %index of the last dro encountered
x_0 = DRO.x_0(k);
Cj = DRO.C_j(k);
Ax = cr3bp.L * (1-cr3bp.mu - x_0);
orbit_2 = init_orbit(cr3bp, cr3bp.l2, cst.orbit.type.PLYAP, cst.orbit.family.PLANAR, Ax, cst);
v_y = sqrt(abs( x_0^2  - Cj + 2*( (1-cr3bp.mu)/abs(x_0+cr3bp.mu) + cr3bp.mu/abs(x_0-1+cr3bp.mu) ) ));
yvg = [ x_0; 0; 0; ;0 ; v_y ; 0];
orbit_2 = orbit_refinement(cr3bp, orbit_2, default, yvg, cst);


%% Continuation


% Waitbar
h = waitbar(0,'First Loop in progress...');

% Loop using the 16 known DRO
MaxIter = 16;
i = length(Dro_data_2.x); % begin after the smaller DROs
while k < MaxIter
    while orbit_2.y0(1) > DRO.x_0(k+1) && k<16
        
        if (orbit_1.y0(1)-orbit_2.y0(1))*cr3bp.L < 100 % to prevent the step decreasing under 100km
            x_0 = Dro_data.x(i-2);
            Cj = Dro_data.C(i-2);
            Ax = cr3bp.L * (1-cr3bp.mu - x_0);
            orbit_1 = init_orbit(cr3bp, cr3bp.l2, cst.orbit.type.PLYAP, cst.orbit.family.PLANAR, Ax, cst);
            v_y = sqrt(abs( x_0^2  - Cj + 2*( (1-cr3bp.mu)/abs(x_0+cr3bp.mu) + cr3bp.mu/abs(x_0-1+cr3bp.mu) ) ));
            yvg = [ x_0; 0; 0; ;0 ; v_y ; 0];
            orbit_1 = orbit_refinement(cr3bp, orbit_1, default, yvg, cst);
        end
        
        yv = orbit_2.y0 + (orbit_2.y0 - orbit_1.y0);
        orbit_1 = orbit_2;
        orbit_2 = orbit_refinement(cr3bp, orbit_2, default, yv, cst);
        Dro_data.T(i) = orbit_2.T;
        x = orbit_2.y0(1);
        Dro_data.x(i) = x;
        Dro_data.C(i) = x^2 + 2*( (1-cr3bp.mu)/abs(x+cr3bp.mu) + cr3bp.mu/abs(x-1+cr3bp.mu) ) - orbit_2.y0(5)^2;
        i = i+1;
    end
    k = k+1;
    waitbar(k/MaxIter);
    x_0 = DRO.x_0(k);
    Cj = DRO.C_j(k);
    Ax = cr3bp.L * (1-cr3bp.mu - x_0);
    orbit_1 = init_orbit(cr3bp, cr3bp.l2, cst.orbit.type.PLYAP, cst.orbit.family.PLANAR, Ax, cst);
    v_y = sqrt(abs( x_0^2  - Cj + 2*( (1-cr3bp.mu)/abs(x_0+cr3bp.mu) + cr3bp.mu/abs(x_0-1+cr3bp.mu) ) ));
    yvg = [ x_0; 0; 0; ;0 ; v_y ; 0];
    orbit_1 = orbit_refinement(cr3bp, orbit_1, default, yvg, cst);
    
end
close(h)

% Loop to continue after the first 16 DRO
% Waitbar
h = waitbar(0,'Second Loop in progress...');
Earth_limit = 30000/cr3bp.L; % not to close from Earth
while orbit_2.y0(1) > Earth_limit
    
    if (orbit_1.y0(1)-orbit_2.y0(1))*cr3bp.L < 100 % to prevent the step decreasing under 100km
        x_0 = Dro_data.x(i-2);
        Cj = Dro_data.C(i-2);
        Ax = cr3bp.L * (1-cr3bp.mu - x_0);
        orbit_1 = init_orbit(cr3bp, cr3bp.l2, cst.orbit.type.PLYAP, cst.orbit.family.PLANAR, Ax, cst);
        v_y = sqrt(abs( x_0^2  - Cj + 2*( (1-cr3bp.mu)/abs(x_0+cr3bp.mu) + cr3bp.mu/abs(x_0-1+cr3bp.mu) ) ));
        yvg = [ x_0; 0; 0; ;0 ; v_y ; 0];
        orbit_1 = orbit_refinement(cr3bp, orbit_1, default, yvg, cst);
    end
    
    yv = orbit_2.y0 + (orbit_2.y0 - orbit_1.y0);
    orbit_1 = orbit_2;
    orbit_2 = orbit_refinement(cr3bp, orbit_2, default, yv, cst);
    Dro_data.T(i) = orbit_2.T;
    x = orbit_2.y0(1);
    Dro_data.x(i) = x;
    Dro_data.C(i) = x^2 + 2*( (1-cr3bp.mu)/abs(x+cr3bp.mu) + cr3bp.mu/abs(x-1+cr3bp.mu) ) - orbit_2.y0(5)^2;
    i = i+1;
    waitbar((DRO.x_0(MaxIter)-orbit_2.y0(1))/(DRO.x_0(MaxIter)-Earth_limit));
end

close(h)
save Dro_data Dro_data


%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
   figure(4);
   view([-47 28]);
end

%%
figure;
hold on
grid on
plot(Dro_data.x, Dro_data.C);

