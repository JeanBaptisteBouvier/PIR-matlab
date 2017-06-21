%--------------------------------------------------------------------------
% Plots an NRO orbit.
%--------------------------------------------------------------------------
%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Selection of LVLH reference primary
% If primary = 1 it will be chosen the ECI (the primary is the Earth)
% otherwise LCI (the primary is the Moon) 

primary = 0;

if(primary == 1)  % It will be chosen the ECI 
    rsyn2rlvlh_c        = @rsyn2rlvlh_e;
else              % It will be chosen the LCI 
    rsyn2rlvlh_c        = @rsyn2rlvlh;
end

%% Selection of Rendezvous strategy 
% If strategy = 1 it will be chosen the Glide approach
% otherwise the Corridor approach

strategy = 0;
if(strategy == 1)  % It will be chosen the Glide approach
    line_of_sight_c        = @line_of_sight_glide;
else              % It will be chosen the Corridor approach
    line_of_sight_c        = @line_of_sight_corridor;
end


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
theta_T = 180*pi/180;
theta_C = 178*pi/180;

% selection of method
choice = [0 0 1 0]; %[1 0 0 0] CW, [0 1 0 0] SL, [0 0 1 0] LR, [1 1 1 0] previous 3, [0 0 0 1] continuation


t_T = interp1(nro_T.alpha, nro_T.tv, theta_T, 'spline');
t_C = interp1(nro_C.alpha, nro_C.tv, theta_C, 'spline');
x0_T   = state_time(t_T,nro_T)';  %initial position of the target in the synodic frame
x0_C   = state_time(t_C,nro_C)';  %initial position of the chaser in the synodic frame
x0_rel = x0_C - x0_T;             %synodic frame 


% Inertial reference frame and comment the line after 
x0_LVLH     = rsyn2rlvlh_c(t_T, x0_rel, x0_T, cr3bp.mu);   %relative vector in LVLH 
x0 =  x0_LVLH(1:3);

% Parameters' setting
alpha = 13*pi/180;
beta = 13*pi/180;
phi = 5*pi/180;
min_dist = 500e-3/cr3bp.L;
TOF = 3600*10;

% selection of hold points and TOFs
out = line_of_sight_c(alpha, beta, phi, TOF, x0, min_dist, cr3bp);
hold_points_dim = out.hold_points; % hold points in the LVLH framework
TOF_dim  = [out.delta_T; 3600];
TOF_cont = TOF_dim./2;
int_time = 500;
dock        = [0 0 0];
r_f = [x0_LVLH(1:3)';hold_points_dim; dock];
np  = size(r_f,1);


% other settings
intervals = 100; % n of intervals in which the TOF is diveded for LR and CW
toll = 1e-12;
it_max = 100;
r_sphere = 2/cr3bp.L;   % radius of the Approach Sphere to ensure trajectory safety policies [km]
r_sphere_KOS = 200e-3/cr3bp.L;   % radius of the Keep Out Sphere to ensure trajectory safety policies [km]
options = odeset('Reltol', default.ode113.RelTol, 'Abstol', default.ode113.AbsTol ,'Events',@(t,y)odezero_sphere_check(t,y,r_sphere)); 
options2 = odeset('Reltol', default.ode113.RelTol, 'Abstol', default.ode113.AbsTol ,'Events',@(t,y)odezero_sphere_check(t,y,r_sphere_KOS)); 

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
settings.options2 = options2;
settings.cr3bp = cr3bp;

% Selection of the hold point for the dispersion. 
%
% If disp_HP = 1, the dispersion will be added to the velocity of the
% initial point that is the chaser iniial velocity;
% if disp_HP = 2 the dispersion will be added to the first HP
% if disp_HP = 0 the dispersion will not be added anywhere
% disp_HP can be chosen between 1 and np-1
% MC_runs represents the number of Monte Carlo simulations that the user
% want to run. MC_runs can be chosen between 0 and inf. 
% ATTENTION : MC_runs = 0 just if disp_HP = 0!

disp_HP = 0;
MC_runs = 0;

if disp_HP > np-1 || disp_HP < 0
        error('The hold point chosen is outside the validity limit')
elseif disp_HP ~= 0 && MC_runs == 0 
        error('The number of Monte Carlo simulations chosen is not valid')
elseif disp_HP == 0 && MC_runs ~= 0
        error('It is not possible to set the number of Monte Carlo simulations different from zero if there is not dispersion')
end

% Inertial reference frame and comment the line after
output_ch = Rendezvous(choice, hold_points_dim, TOF_dim, init_cond, orbits, settings, continuation,1,0, primary, disp_HP, MC_runs);




