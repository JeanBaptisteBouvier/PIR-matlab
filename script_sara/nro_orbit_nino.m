%--------------------------------------------------------------------------
% Plots an NRO orbit.
%
% BLB 2016
%--------------------------------------------------------------------------
%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Change of parameters wrt default
default.plot.XY          = false;  %plot also the results in X-Y plane
default.plot.XZ          = false;  %plot also the results in X-Z plane
default.plot.YZ          = false;  %plot also the results in Y-Z plane
default.plot.TD          = false;   %plot also the results in three-D
%% Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%% Orbit TARGET
%Initialize NRO 
nro_T = init_nro(cr3bp, cr3bp.l2, cst.orbit.family.SOUTHERN, cst);

% Interpolation and plot
nro_T = nro_interpolation(cr3bp, nro_T, nro_init_EML2, default, cst, 'Az', 70000); %Az o periselene??

%% Orbit CHASER
%Initialize NRO 
nro_C = init_nro(cr3bp, cr3bp.l2, cst.orbit.family.SOUTHERN, cst);

% Interpolation and plot
nro_C = nro_interpolation(cr3bp, nro_C, nro_init_EML2, default, cst, 'Az', 69960);


%% Nonlinear Relative dynamics in synodic frame

% initial conditions
% theta_T = 350*pi/180;
% theta_C = 315*pi/180;

theta_T = 180*pi/180;
theta_C = 178*pi/180;

% selection of method
choice = [1 0 0 0]; %[1 0 0 0] CW, [0 1 0 0] SL, [0 0 1 0] LR, [1 1 1 0] previous 3, [0 0 0 1] continuation

% selection of hold points and TOFs
% hold_points_dim = [-5 0 -1];
% TOF_dim  = [3600*2 3600];
% TOF_cont = [3600 600]; 

t_T = interp1(nro_T.alpha, nro_T.tv, theta_T, 'spline');
t_C = interp1(nro_C.alpha, nro_C.tv, theta_C, 'spline');
x0_T   = state_time(t_T,nro_T)';  %initial position of the target in the synodic frame
x0_C   = state_time(t_C,nro_C)';  %initial position of the chaser in the synodic frame
x0_rel = x0_C - x0_T;             %synodic frame 
x0_LVLH     = rsyn2rlvlh(t_T, x0_rel, x0_T, cr3bp.mu);   %relative vector in LVLH with z-axis toward the Moon.
x0 = x0_LVLH(1:3)*cr3bp.L;

alpha = 35*pi/180;
beta = 35*pi/180;
phi = 0.5*pi/180;
min_dist = 0.5;
TOF = 3600*10;
% out =  line_of_sight_corridor( alpha, beta, phi, TOF , x0, min_dist);
out = line_of_sight_glide(alpha, phi, TOF, x0, min_dist);
hold_points_dim = out.hold_points;
% hold_points_dim = [-5 1 0]; % control point in the LVLH framework
TOF_dim  = [out.delta_T; 3600];
TOF_cont = TOF_dim./2;
int_time = 300;


% other settings
intervals = 100; %n di intervalli in cui è diviso il TOF per Luquette & CW
toll = 1e-12;
it_max = 100;
options = odeset('Reltol', default.ode113.RelTol, 'Abstol', default.ode113.AbsTol);

%structures as inputs
init_cond.T = theta_T;
init_cond.C = theta_C;
continuation.TOF_cont = TOF_cont;
continuation.int = int_time;
orbits.T = nro_T;
orbits.C = nro_C;

settings.Luq_int = intervals;
settings.toll = toll;
settings.it_max = it_max;
settings.options = options;
settings.cr3bp = cr3bp;

output_ch = Rendezvous_choice_l(choice, hold_points_dim, TOF_dim, init_cond, orbits, settings, continuation, 1);




