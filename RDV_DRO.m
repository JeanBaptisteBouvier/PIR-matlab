% Rdv between two DROs
% 
% Origin position on a DRO : theta = 0, on X-axis, closest from the Earth
% Theta is growing clockwise
% 
% 
% 
% By JBB

init;

%% Orbits of the Target and the Chaser in the synodic frame

% Orbit of the Target
[dro_T, default, ~] = A_DRO(99900, default, cst);
theta_T = 140*pi/180;                   % initial positions of the Target


% Orbit of the Chaser
[dro_C, default, cr3bp] = A_DRO(100000, default, cst);
theta_C = 130*pi/180;                   % initial positions of the Chaser


% plot
% orbit_plot(dro_T, default, 'red')
% orbit_plot(dro_C, default, 'blue')

%% Initialisation of the main structures

% selection of hold points and TOFs
% if no hold point is required, then the following vector can be set to []
t_T = interp1(dro_T.theta, dro_T.tv, theta_T, 'spline');
t_C = interp1(dro_C.theta, dro_C.tv, theta_C, 'spline');
x0_T   = state_time(t_T,dro_T)';  %initial position of the target in the synodic frame
x0_C   = state_time(t_C,dro_C)';  %initial position of the chaser in the synodic frame
x0_rel = x0_C - x0_T;             %synodic frame 
x0_LVLH     = rsyn2rlvlh(t_T, x0_rel, x0_T, cr3bp.mu);   %relative vector in LVLH with z-axis toward the Moon.
x0 = x0_LVLH(1:3)*cr3bp.L;

alpha = 10*pi/180;
beta = 10*pi/180;
phi = 7.5*pi/180;
theta = 0;
min_dist = 0.5;
TOF = 3600*10;
out = line_of_sight_glide(alpha, phi, theta, TOF, x0, min_dist);
hold_points_dim = out.hold_points;
% hold_points_dim = [-5 1 0]; % control point in the LVLH framework
TOF_dim  = [out.delta_T; 3600];
TOF_cont = TOF_dim./2; 
int_time = 300;

% other settings
intervals = 100; % number of time intervals in TOF 
it_max = 100; % iteration max
options = odeset('Reltol', default.ode113.RelTol, 'Abstol', default.ode113.AbsTol);

%structures as inputs
init_cond.T = theta_T;
init_cond.C = theta_C;
continuation.TOF_cont = TOF_cont;
continuation.int = int_time;
orbits.T = dro_T;
orbits.C = dro_C;
settings.Luq_int = intervals;
settings.it_max = it_max;
settings.options = options;
settings.cr3bp = cr3bp;
settings.toll = 1e-12;

%% Computation of the rendezvous trajectory for a user-given number of hold points, initial orbits and TOFS
output = Rendezvous_choice_dro(hold_points_dim, TOF_dim, init_cond, orbits, settings, continuation);

sol = output.linear;
Lam = output.Lambert;
dV = output.dV;

