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
nro_C = nro_interpolation(cr3bp, nro_C, nro_init_EML2, default, cst, 'Az', 69900);


%% Nonlinear Relative dynamics in synodic frame

% initial conditions

theta_T = [0 180 90]*pi/180;
step = 0.5*pi/180;
d_theta = [ 10 0.33 3.57]*pi/180;
theta_C_p1 = (theta_T(1)- d_theta(1) + 360*pi/180): step : 360*pi/180;
theta_C_p2 = 0.5*pi/180 :step: d_theta(1);
theta_C_p = [ theta_C_p1  theta_C_p2];
theta_C_a = (theta_T(2) - d_theta(2)):step:(theta_T(2) + d_theta(2));
theta_C_i = (theta_T(3) - d_theta(3)):step:(theta_T(3) + d_theta(3));
theta_C = [theta_C_p  theta_C_a  theta_C_i];
TOF = [600:300:5700];

% selection of method
choice = [1 1 1 0]; %[1 0 0 0] CW, [0 1 0 0] SL, [0 0 1 0] LR, [1 1 1 0] previous 3, [0 0 0 1] continuation

% selection of hold points and TOFs
hold_points_dim = [];
TOF_cont = 3600*2; 
int_time = 300;

% other settings
intervals = 100; %n di intervalli in cui è diviso il TOF per Luquette & CW
toll = 1e-12;
it_max = 100;
options = odeset('Reltol', default.ode113.RelTol, 'Abstol', default.ode113.AbsTol);

% variables initialization 
n = length(theta_C);
m = length(TOF);
p = length(theta_T);

err_CW_e = zeros(n,m);
err_SL_e = zeros(n,m);
err_LR_e = zeros(n,m);
x_init = zeros(n,m);
y_init = zeros(n,m);
z_init = zeros(n,m);
init_dist = zeros(n,m);
err_CW_P_e = zeros(length(theta_C_p) , m);
err_SL_P_e = zeros(length(theta_C_p) , m);
err_LR_P_e = zeros(length(theta_C_p) , m);

err_CW_A_e = zeros(length(theta_C_a) , m);
err_SL_A_e = zeros(length(theta_C_a) , m);
err_LR_A_e = zeros(length(theta_C_a) , m);

err_CW_I_e = zeros(length(theta_C_i) , m);
err_SL_I_e = zeros(length(theta_C_i) , m);
err_LR_I_e = zeros(length(theta_C_i) , m);

% Errors compitation
for k = 1:p
  for i = 1:n
    for j = 1:m
    %structures as inputs
    init_cond.T = theta_T(k);
    init_cond.C = theta_C(i);
    continuation.TOF_cont = TOF_cont;
    continuation.int = int_time;
    orbits.T = nro_T;
    orbits.C = nro_C;
    TOF_dim  = TOF(j);

    settings.Luq_int = intervals;
    settings.toll = toll;
    settings.it_max = it_max;
    settings.options = options;
    settings.cr3bp = cr3bp;

    output_ch = Rendezvous_errors_e(choice, hold_points_dim, TOF_dim, init_cond, orbits, settings, continuation, 0);
     
    switch k
        case 1 % Periselene
            err_CW_P_e(i,j)    = output_ch.linear{1,1}.CW.err1.*cr3bp.L;
            err_SL_P_e(i,j)    = output_ch.linear{1,1}.SL.err1.*cr3bp.L;
            err_LR_P_e(i,j)    = output_ch.linear{1,1}.LR.err1.*cr3bp.L;
            x_init_P_e(i,j)    = output_ch.init_LVLH(1)*cr3bp.L;
            y_init_P_e(i,j)    = output_ch.init_LVLH(2)*cr3bp.L;
            z_init_P_e(i,j)    = output_ch.init_LVLH(3)*cr3bp.L;
            init_pos_P_e       = [x_init_P_e(i,j) y_init_P_e(i,j) z_init_P_e(i,j)];
            init_dist_P_e(i,j) = norm(init_pos_P_e);
        
        case 2 % Aposelene
            err_CW_A_e(i,j)    = output_ch.linear{1,1}.CW.err1.*cr3bp.L;
            err_SL_A_e(i,j)    = output_ch.linear{1,1}.SL.err1.*cr3bp.L;
            err_LR_A_e(i,j)    = output_ch.linear{1,1}.LR.err1.*cr3bp.L;
            x_init_A_e(i,j)    = output_ch.init_LVLH(1)*cr3bp.L;
            y_init_A_e(i,j)    = output_ch.init_LVLH(2)*cr3bp.L;
            z_init_A_e(i,j)    = output_ch.init_LVLH(3)*cr3bp.L;
            init_pos_A_e       = [x_init_A_e(i,j) y_init_A_e(i,j) z_init_A_e(i,j)];
            init_dist_A(i,j) = norm(init_pos_A_e);
            
        case 3  % Intermediate 
            err_CW_I_e(i,j)    = output_ch.linear{1,1}.CW.err1.*cr3bp.L;
            err_SL_I_e(i,j)    = output_ch.linear{1,1}.SL.err1.*cr3bp.L;
            err_LR_I_e(i,j)    = output_ch.linear{1,1}.LR.err1.*cr3bp.L;
            x_init_I_e(i,j)    = output_ch.init_LVLH(1)*cr3bp.L;
            y_init_I_e(i,j)    = output_ch.init_LVLH(2)*cr3bp.L;
            z_init_I_e(i,j)    = output_ch.init_LVLH(3)*cr3bp.L;
            init_pos_I_e       = [x_init_I_e(i,j) y_init_I_e(i,j) z_init_I_e(i,j)];
            init_dist_I(i,j) = norm(init_pos_I_e); 
    end
    end
  end
end
err_CW_e = [err_CW_P_e; err_CW_A_e; err_CW_I_e];
err_SL_e = [err_SL_P_e; err_SL_A_e; err_SL_I_e];
err_LR_e = [err_LR_P_e; err_LR_A_e; err_LR_I_e];

%% Plot

% MAX ERROR
figure(1)
x1 = [1:3];
x2 = [30:32];
x3 = [57:59];
x = [ x1; x2; x3];
y_CW_MAX_P_e = max(max(err_CW_P_e));
y_SL_MAX_P_e = max(max(err_SL_P_e));
y_LR_MAX_P_e = max(max(err_LR_P_e));
y_MAX_P_e = [y_CW_MAX_P_e y_SL_MAX_P_e  y_LR_MAX_P_e];

y_CW_MAX_A_e = max(max(err_CW_A_e));
y_SL_MAX_A_e = max(max(err_SL_A_e));
y_LR_MAX_A_e = max(max(err_LR_A_e));
y_MAX_A_e = [y_CW_MAX_A_e y_SL_MAX_A_e  y_LR_MAX_A_e];


y_CW_MAX_I_e = max(max(err_CW_I_e));
y_SL_MAX_I_e = max(max(err_SL_I_e));
y_LR_MAX_I_e = max(max(err_LR_I_e));
y_MAX_I_e = [y_CW_MAX_I_e y_SL_MAX_I_e  y_LR_MAX_I_e];

y_max_e = [y_MAX_P_e ;y_MAX_A_e; y_MAX_I_e];  

bar(x ,y_max_e)
title('MAX Dimensional Error [km]')
legend('CW', 'SL' , 'LR')
name = { 'Periselene' , 'Aposelene' , 'Intermediate'};
set(gca,'xticklabel',name)

% MIN ERROR
figure(2)
y_CW_MIN_P_e = min(min(err_CW_P_e));
y_SL_MIN_P_e = min(min(err_SL_P_e));
y_LR_MIN_P_e = min(min(err_LR_P_e));
y_MIN_P_e = [y_CW_MIN_P_e y_SL_MIN_P_e  y_LR_MIN_P_e];

y_CW_MIN_A_e = min(min(err_CW_A_e));
y_SL_MIN_A_e = min(min(err_SL_A_e));
y_LR_MIN_A_e = min(min(err_LR_A_e));
y_MIN_A_e = [y_CW_MIN_A_e y_SL_MIN_A_e  y_LR_MIN_A_e];


y_CW_MIN_I_e = min(min(err_CW_I_e));
y_SL_MIN_I_e = min(min(err_SL_I_e));
y_LR_MIN_I_e = min(min(err_LR_I_e));
y_MIN_I_e = [y_CW_MIN_I_e y_SL_MIN_I_e  y_LR_MIN_I_e];

y_min_e = [y_MIN_P_e; y_MIN_A_e; y_MIN_I_e];
bar(x,y_min_e)
title('MIN Dimensional Error [km]')
legend('CW', 'SL' , 'LR')
set(gca,'xticklabel',name)

% MEAN ERROR
figure(3)
y_CW_MEAN_P_e = mean2(err_CW_P_e);
y_SL_MEAN_P_e = mean2(err_SL_P_e);
y_LR_MEAN_P_e = mean2(err_LR_P_e);
y_MEAN_P_e = [y_CW_MEAN_P_e y_SL_MEAN_P_e  y_LR_MEAN_P_e];

y_CW_MEAN_A_e = mean2(err_CW_A_e);
y_SL_MEAN_A_e = mean2(err_SL_A_e);
y_LR_MEAN_A_e = mean2(err_LR_A_e);
y_MEAN_A_e = [y_CW_MEAN_A_e y_SL_MEAN_A_e  y_LR_MEAN_A_e];


y_CW_MEAN_I_e = mean2(err_CW_I_e);
y_SL_MEAN_I_e = mean2(err_SL_I_e);
y_LR_MEAN_I_e = mean2(err_LR_I_e);
y_MEAN_I_e = [y_CW_MEAN_I_e y_SL_MEAN_I_e  y_LR_MEAN_I_e];

y_mean_e = [ y_MEAN_P_e; y_MEAN_A_e; y_MEAN_I_e];
bar(x,y_mean_e)
title('MEAN Dimensional Error [km]')
legend('CW', 'SL' , 'LR')
set(gca,'xticklabel',name)


