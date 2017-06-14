function output = Rendezvous_choice_dro(hold_points_dim, TOF_dim, init_cond, orbits, settings, continuation)
% Author: Antonino Campolo 2016/2017 (with the Holy help of the Fabolous Mr. BLB) 

% This code computes the rendezvous trajectory for a user-given number of
% hold points, initial orbits and TOFs 

% hold_points_dim is a matrix whose rows are the LVLH DIMENSIONAL (km) positions of the hold
% points. If no hold points are needed, hold_points_dim =[];

% TOF_dim is a vector which contains the (n. of hold points +1) DIMENSIONAL (s) time of flights for
% each transfer from one hold point to the next. At least one TOF is needed
% to compute the rendezvous

% theta_T & theta_C (radians) are the "true anomalies" of the target and the chaser
% on the planar orbit resulting from the projection of the 3D chosen orbit
% on the plane formed by its position and velocity vector at the periselene

% dro_T & dro_C are the structures regarding target and chaser's orbits
% (from main)

% toll is the tolerance on the error in the Lambert cycle. The error is
% computed as the norm of the vectorial difference between the integrated
% final position and the desired final position

% it_max is the maximum number of iteration allowed for the Lambert's cycle

% intervals is the number of intervals in which the Luquette's transfer is
% divided. For each interval, the STM is evaluated in order to compose a
% final STM which is better conditioned

% TOF_cont is a vector which contains the starting DIMENSIONAL (s) Time of 
% Flights for each transfer for the continuation method. 

% int_time is the DIMENSIONAL (s) time separation between two evaluation of
% Lambert's arcs in the continuation algorithm
% choice is a 3x1 vector wich selects a particular first guess for the


% choice = [1 0 0 0], CW is selected.
% choice = [0 1 0 0], SL is selected
% choice = [0 0 1 0], LR is selected
% choice = [1 1 1 0], all three previous methods are used to generate an
% hybrid solution
% choice = [0 0 0 1], continuation algorithm is selected
% No other combinations are possible

choice = [1 1 1 0]; 

%% Inputs generation

theta_T = init_cond.T;
theta_C = init_cond.C;
dro_T = orbits.T;
dro_C = orbits.C;
cr3bp = settings.cr3bp;

TOF_cont = continuation.TOF_cont;
int_time = continuation.int;

t_T = interp1(dro_T.theta, dro_T.tv, theta_T, 'spline'); %defines the initial relative position of synodic/inertial frame
t_C = interp1(dro_C.theta, dro_C.tv, theta_C, 'spline');
                   
x0_T   = state_time(t_T,dro_T)';  %initial position of the target in the synodic frame
x0_C   = state_time(t_C,dro_C)';  %initial position of the chaser in the synodic frame
x0_rel = x0_C - x0_T;             %synodic frame 
phi0   = reshape(eye(6,6),[],1);
 
% states in LVLH
x0_LVLH     = rsyn2rlvlh(t_T, x0_rel, x0_T, cr3bp.mu);   %relative vector in LVLH with z-axis toward the Moon.
hold_points = hold_points_dim./cr3bp.L;                %adimensional distance for hold points in LVLH
dock        = [0 0 0];                                 %position of the target in LVLH frame


% successive positions and velocities in the LVLH frame
r_f = [x0_LVLH(1:3)'; hold_points; dock];  
np  = size(r_f,1); % nbr of hold points + initial point + final point
v_f = [x0_LVLH(4:6)'; zeros(np-1,3)];  
 
% TOF and times vector
TOF = TOF_dim.*2*pi/cr3bp.T;

% possible errors
if length(TOF) ~= np-1
   error('The number of TOFs must be equal to the number of hold points + 1')
end


% build a time interval containing the discrete time inervals station with
% the initial shift of the synodic frame
t_f    = zeros(np,1);
t_f(1) = t_T;
for i=1:np-1
    t_f(i+1) = t_f(i) + TOF(i);  %t_f must be adimensional
end

% states in Synodic frame
m_hp   = size(hold_points,1);
hp_syn = zeros(m_hp,6);

if m_hp > 0
   for i=1:m_hp
       % CoC of the hold_points from the LVLH into the Synodic frame
       hp_syn(i,:) = rlvlh2rsyn(t_f(i+1), [hold_points(i,:) 0 0 0]', state_time(t_f(i+1),dro_T)', cr3bp.mu);
   end
   
   r_f_syn = [x0_rel(1:3)'; hp_syn(:,1:3); dock];
   v_f_syn = [x0_rel(4:6)'; zeros(np-1,3)];
 
else
   hp_syn = []; 
   r_f_syn = [x0_rel(1:3)'; hp_syn; dock];
   v_f_syn = [x0_rel(4:6)'; zeros(np-1,3)]; 
end
 

%% First guesses & Lambert

% np = 2 when no hold points
sol = cell(1,np-1);
for i = 1:np-1
   
    LVLH.r         = r_f(i:i+1,:);
    LVLH.v         = v_f(i,:);
    SYN.r          = r_f_syn(i:i+1,:); 
    SYN.v          = v_f_syn(i,:);
    SYN.r_next     = SYN.r(2,:);
    times_lin.T    = TOF(i); 
    times_lin.tt   = t_f(i);
    times_lin.tt_2 = t_f(i+1);
    times_cont.TOF = TOF_cont(i);
    times_cont.int = int_time;
    
    pop = Targeting(choice, LVLH, SYN, times_lin, times_cont, settings, dro_T, phi0, i);
    sol{i} = pop.sol;
end


 %% Building Hybrid Solution
% Lam = building(choice, sol, np, it_max);
 
%% deltaV calculation
Lam = cell(1,np-1);
for i=1:np-1
    Lam{i} = sol{i}.Lam;
end
dV = dV_eval(Lam, x0_rel, np, cr3bp);

%% Outputs
% solution
output.linear = sol;
% deltaV
output.dV = dV;
 
%% Plotting

% Choice is the same vector for each hold points
plot_rdv_DRO(choice, sol, r_f_syn, r_f, np, cr3bp, dro_T, dro_C, x0_T, x0_C, t_f, t_C);




end