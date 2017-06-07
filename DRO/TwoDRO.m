% Rdv between two DROs
% 
% Origin position on a DRO : theta = 0°, on X-axis, closest from the Earth
% Theta is growing clockwise
% 
% 
% 
% By JBB

init;

%% Orbits

% Orbit of the Target
[dro_T, default, ~] = A_DRO(40000, default, cst);
theta_T = 180*pi/180;                   % initial positions of the Target


% Orbit of the Chaser
[dro_C, default, cr3bp] = A_DRO(100000, default, cst);
theta_C = 178*pi/180;                   % initial positions of the Chaser


% plot
% orbit_plot(dro_T, default, 'red')
% orbit_plot(dro_C, default, 'blue')

%% 

% selection of hold points and TOFs
hold_points_dim = [-5 0 -1]; % ???
TOF_dim  = [3600*2 3600];
TOF_cont = [3600 600]; 
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

% output_ch = Rendezvous_choice2(hold_points_dim, TOF_dim, init_cond, orbits, settings, continuation);


