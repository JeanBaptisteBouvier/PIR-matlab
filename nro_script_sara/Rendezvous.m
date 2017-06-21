function output = Rendezvous (choice, hold_points_dim, TOF_dim, init_cond, orbits, settings, continuation, plot_yes, plot_fg, primary, disp_HP, MC_runs)
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

% nro_T & nro_C are the structures regarding target and chaser's orbits
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
% Lambert's cycle. If 
% choice = [1 0 0 0], CW is selected.
% choice = [0 1 0 0], SL is selected
% choice = [0 0 1 0], LR is selected
% choice = [1 1 1 0], all three previous methos are used to generate an
% hybrid solution
% choice = [0 0 0 1], continuation algorithm is selected
% No other combinations are possible

% if plot is required, plot_yes = 1. Otherwise, put any other number

%% Selection of LVLH reference primary

if(primary == 1)  % It will be choosen the ECI 
    rsyn2rlvlh_c        = @rsyn2rlvlh_e;
    rlvlh2rsyn_c        = @rlvlh2rsyn_e;
else              % It will be choosen the LCI 
    rsyn2rlvlh_c        = @rsyn2rlvlh;
    rlvlh2rsyn_c        = @rlvlh2rsyn;
end


%% Inputs generation

theta_T = init_cond.T;
theta_C = init_cond.C;
TOF_cont = continuation.TOF_cont;
int_time = continuation.int;
nro_T = orbits.T;
nro_C = orbits.C;
cr3bp = settings.cr3bp;

t_T = interp1(nro_T.alpha,nro_T.tv,theta_T,'spline');
t_C = interp1(nro_C.alpha,nro_C.tv,theta_C,'spline');

t      = t_T;                     %defines the initial relative position of synodic/inertial frame
x0_T   = state_time(t_T,nro_T)';
x0_C   = state_time(t_C,nro_C)';
x0_rel = x0_C - x0_T;             %synodic frame 
phi0   = reshape(eye(6,6),[],1);
 
% states in LVLH
x0_LVLH = rsyn2rlvlh_c(t, x0_rel, x0_T, cr3bp.mu);   %adim
hold_points = hold_points_dim; %              %LVLH
dock        = [0 0 0];

% successive positions and velocities in the LVLH frame
r_f = [x0_LVLH(1:3)';hold_points_dim; dock];
np  = size(r_f,1);

v_f = [x0_LVLH(4:6)'; zeros(np-1,3)];  
 
% TOF and times vector
TOF = TOF_dim.*2*pi/cr3bp.T;


% possible errors
if length(TOF) ~= np-1
   error('The number of TOFs must be equal to the n.of hold points + 1')
end

if length(TOF_cont) ~= np-1
   error('The number of TOFs for the continuation method must be equal to the n.of hold points + 1')
end

for i=1:np-1
    if TOF_cont(i).*2*pi/cr3bp.T >= TOF(i)
       error('TOF_cont must be < than TOF for each transfer')
    end
end

% times for each hold point
t_f    = zeros(np,1);
t_f(1) = t;   
for i=1:np-1
    t_f(i+1) = t+sum(TOF(1:i));  %t must be adimensional
end

% states in Synodic frame
m_hp   = size(hold_points,1);
hp_syn = zeros(m_hp,6);

if m_hp > 0
   for i=1:m_hp
       hp_syn(i,:) = rlvlh2rsyn_c(t_f(i+1), [hold_points(i,:) 0 0 0]', state_time(t_f(i+1),nro_T)', cr3bp.mu);
   end
   
   r_f_syn = [x0_rel(1:3)'; hp_syn(:,1:3); dock];
   v_f_syn = [x0_rel(4:6)'; zeros(np-1,3)];
 
else
   hp_syn = []; 
   r_f_syn = [x0_rel(1:3)'; hp_syn; dock];
   v_f_syn = [x0_rel(4:6)'; zeros(np-1,3)]; 
end
 
%% First guesses & Lambert
if disp_HP == 0
    sol = cell(1,np-1);
    Lam = cell(1,np-1);
else
    Lam_t_HP = cell(1,MC_runs);
    Lam_HP = cell(1,MC_runs);
end

it = 1;
t_events_AS = cell(1,MC_runs+1);
y_events_AS = cell(1,MC_runs+1);
t_events_KOS = cell(1,MC_runs+1);
y_events_KOS = cell(1,MC_runs+1);
sol_HP = cell(1,disp_HP+MC_runs);
for i = 1:np-1
    r_LVLH = r_f(i:i+1,:);
    v_LVLH = v_f(i,:);
    r_syn  = r_f_syn(i:i+1,:);
    r_next = r_syn(2,:);
    v_syn  = v_f_syn(i,:);
    TOF_ca = TOF_cont(i);
    T      = TOF(i);
    tt     = t_f(i);
    tt_2   = t_f(i+1);
    
    LVLH.r         = r_LVLH;
    LVLH.v         = v_LVLH;
    SYN.r          = r_syn; 
    SYN.v          = v_syn;
    SYN.r_next     = r_next;
    times_lin.T    = T; 
    times_lin.tt   = tt;
    times_lin.tt_2 = tt_2;
    times_cont.TOF = TOF_ca;
    times_cont.int = int_time;
    n_transf = i;

    pop = Terminator(choice, LVLH, SYN, times_lin, times_cont, settings, nro_T, phi0, n_transf, primary);
    sol{i} = pop.sol;
    
    %    % Building Lambert
     Lam{i} = sol{i}.Lam;
     Lam_t{i} = sol{i}.Lam_t;
         
    % Dispersion propagation
    if i == disp_HP
        if disp_HP == 1
            V_minus = x0_rel(4:6)';
        else
            V_minus = Lam{i-1}(end,4:6);
        end
        V_plus = Lam{i}(1,4:6);
        

        while it <= MC_runs
            
            [ Lam_t_HP{it}, Lam_HP{it}, t_events_AS{it}, y_events_AS{it}, t_events_KOS{it}, y_events_KOS{it}] = Dispersion(V_plus, V_minus,  Lam{i}(1,1:3), Lam_t{i}(1), state_time(Lam_t{i}(1), nro_T),  24*3600 , 0.01, 1e-3*180/pi, phi0, cr3bp, settings.options, settings.options2);
             sol_HP{it+i}.Lam_HP = Lam_HP{it};
             sol_HP{it+i}.Lam_t_HP = Lam_t_HP{it};
             it = it+1;
        end
        
      break
    end
end


%% deltaV calculation
dV = dV_comp(Lam, x0_rel, np, cr3bp, disp_HP, MC_runs, sol_HP);

% Propagation of the Lambert arc in the disp_HP for a time of 24h without dispersion
if disp_HP ~= 0
    [time, transfer, t_events_n_AS, y_events_n_AS, ~] = ode113(@(t,y)crtbp(t,y,cr3bp.mu),[Lam_t{disp_HP}(1) 24*3600*2*pi/cr3bp.T+Lam_t{disp_HP}(1)],[Lam{disp_HP}(1,1:6), state_time(Lam_t{disp_HP}(1), nro_T), phi0']', settings.options); 
    if t_events_n_AS ~= 0
       [time, transfer, t_events_n_KOS, y_events_n_KOS, ~] = ode113(@(t,y)crtbp(t,y,cr3bp.mu),[Lam_t{disp_HP}(1) 24*3600*2*pi/cr3bp.T+Lam_t{disp_HP}(1)],[Lam{disp_HP}(1,1:6), state_time(Lam_t{disp_HP}(1), nro_T), phi0']', settings.options2); 
       t_events_KOS{MC_runs+1} = t_events_n_KOS;
       y_events_KOS{MC_runs+1} = y_events_n_KOS;
    end
    t_events_AS{MC_runs+1} = t_events_n_AS;
    y_events_AS{MC_runs+1} = y_events_n_AS;
    clear Lam_t{disp_HP}
    clear Lam{disp_HP}
    Lam_t{disp_HP} = time;
    Lam{disp_HP} =  transfer;
end

%% Outputs
output.linear = sol;
output.Lambert = Lam;
output.t_events_AS = t_events_AS;
output.y_events_AS = y_events_AS;
output.t_events_KOS = t_events_KOS;
output.y_events_KOS = y_events_KOS;

% deltaV
output.dV = dV;
 
%% Plotting

if plot_yes == 1
    if disp_HP ==0
         l = 1:np-1;
    else
        l = 1:disp_HP;
     end
    
    for i = l
      if  i ~= disp_HP
           Lam_plot = Lam{i};
           Lam_plot_t = Lam_t{i};
           plot_choice(Lam_plot, Lam_plot_t, r_f_syn, r_f , np, cr3bp, nro_T, nro_C, x0_T, x0_C,t_f,t_C, primary , disp_HP, 'b', 'RDV Trajectory')  %Lambert
      else
           Lam_plot = Lam{i};
           Lam_plot_t = Lam_t{i};
           plot_choice(Lam_plot, Lam_plot_t, r_f_syn, r_f , np, cr3bp, nro_T, nro_C, x0_T, x0_C, t_f,t_C, primary , disp_HP , 'b', 'RDV Trajectory')  %Lambert
           for s = 1: MC_runs
                Lam_plot = Lam_HP{s};
                Lam_plot_t = Lam_t_HP{s};
                plot_choice(Lam_plot, Lam_plot_t, r_f_syn, r_f , np, cr3bp, nro_T, nro_C, x0_T, x0_C,t_f,t_C, primary , disp_HP, 'r', 'Dispersion Trajectory')  %Lambert
           end
      end
    end
end
%% Plot First Guess

if plot_fg ==1 && disp_HP == 0 
    plot_choice_fg(choice, sol, r_f_syn, r_f, np, cr3bp, nro_T, nro_C, x0_T, x0_C, t_f, t_C,primary)
end
