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

err_CW = zeros(n,m);
err_SL = zeros(n,m);
err_LR = zeros(n,m);
x_init = zeros(n,m);
y_init = zeros(n,m);
z_init = zeros(n,m);
init_dist = zeros(n,m);
err_CW_P = zeros(length(theta_C_p) , m);
err_SL_P = zeros(length(theta_C_p) , m);
err_LR_P = zeros(length(theta_C_p) , m);

err_CW_A = zeros(length(theta_C_a) , m);
err_SL_A = zeros(length(theta_C_a) , m);
err_LR_A = zeros(length(theta_C_a) , m);

err_CW_I = zeros(length(theta_C_i) , m);
err_SL_I = zeros(length(theta_C_i) , m);
err_LR_I = zeros(length(theta_C_i) , m);

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

    output_ch = Rendezvous_errors_sara(choice, hold_points_dim, TOF_dim, init_cond, orbits, settings, continuation, 0);
     
    switch k
        case 1 % Periselene
            err_CW_P(i,j)    = output_ch.linear{1,1}.CW.err1.*cr3bp.L;
            err_SL_P(i,j)    = output_ch.linear{1,1}.SL.err1.*cr3bp.L;
            err_LR_P(i,j)    = output_ch.linear{1,1}.LR.err1.*cr3bp.L;
            x_init_P(i,j)    = output_ch.init_LVLH(1)*cr3bp.L;
            y_init_P(i,j)    = output_ch.init_LVLH(2)*cr3bp.L;
            z_init_P(i,j)    = output_ch.init_LVLH(3)*cr3bp.L;
            init_pos_P       = [x_init_P(i,j) y_init_P(i,j) z_init_P(i,j)];
            init_dist_P(i,j) = norm(init_pos_P);
        
        case 2 % Aposelene
            err_CW_A(i,j)    = output_ch.linear{1,1}.CW.err1.*cr3bp.L;
            err_SL_A(i,j)    = output_ch.linear{1,1}.SL.err1.*cr3bp.L;
            err_LR_A(i,j)    = output_ch.linear{1,1}.LR.err1.*cr3bp.L;
            x_init_A(i,j)    = output_ch.init_LVLH(1)*cr3bp.L;
            y_init_A(i,j)    = output_ch.init_LVLH(2)*cr3bp.L;
            z_init_A(i,j)    = output_ch.init_LVLH(3)*cr3bp.L;
            init_pos_A       = [x_init_A(i,j) y_init_A(i,j) z_init_A(i,j)];
            init_dist_A(i,j) = norm(init_pos_A);
            
        case 3  % Intermediate 
            err_CW_I(i,j)    = output_ch.linear{1,1}.CW.err1.*cr3bp.L;
            err_SL_I(i,j)    = output_ch.linear{1,1}.SL.err1.*cr3bp.L;
            err_LR_I(i,j)    = output_ch.linear{1,1}.LR.err1.*cr3bp.L;
            x_init_I(i,j)    = output_ch.init_LVLH(1)*cr3bp.L;
            y_init_I(i,j)    = output_ch.init_LVLH(2)*cr3bp.L;
            z_init_I(i,j)    = output_ch.init_LVLH(3)*cr3bp.L;
            init_pos_I       = [x_init_I(i,j) y_init_I(i,j) z_init_I(i,j)];
            init_dist_I(i,j) = norm(init_pos_I); 
    end
    end
  end
end
err_CW = [err_CW_P; err_CW_A; err_CW_I];
err_SL = [err_SL_P; err_SL_A; err_SL_I];
err_LR = [err_LR_P; err_LR_A; err_LR_I];

%% Plot

% MAX ERROR
figure(1)
x1 = [1:3];
x2 = [30:32];
x3 = [57:59];
x = [ x1; x2; x3];
y_CW_MAX_P = max(max(err_CW_P));
y_SL_MAX_P = max(max(err_SL_P));
y_LR_MAX_P = max(max(err_LR_P));
y_MAX_P = [y_CW_MAX_P y_SL_MAX_P  y_LR_MAX_P];

y_CW_MAX_A = max(max(err_CW_A));
y_SL_MAX_A = max(max(err_SL_A));
y_LR_MAX_A = max(max(err_LR_A));
y_MAX_A = [y_CW_MAX_A y_SL_MAX_A  y_LR_MAX_A];


y_CW_MAX_I = max(max(err_CW_I));
y_SL_MAX_I = max(max(err_SL_I));
y_LR_MAX_I = max(max(err_LR_I));
y_MAX_I = [y_CW_MAX_I y_SL_MAX_I  y_LR_MAX_I];

y_max = [y_MAX_P ;y_MAX_A; y_MAX_I];  

bar(x ,y_max)
title('MAX Dimensional Error [km]')
legend('CW', 'SL' , 'LR')
name = { 'Periselene' , 'Aposelene' , 'Intermediate'}
set(gca,'xticklabel',name)

% MIN ERROR
figure(2)
y_CW_MIN_P = min(min(err_CW_P));
y_SL_MIN_P = min(min(err_SL_P));
y_LR_MIN_P = min(min(err_LR_P));
y_MIN_P = [y_CW_MIN_P y_SL_MIN_P  y_LR_MIN_P];

y_CW_MIN_A = min(min(err_CW_A));
y_SL_MIN_A = min(min(err_SL_A));
y_LR_MIN_A = min(min(err_LR_A));
y_MIN_A = [y_CW_MIN_A y_SL_MIN_A  y_LR_MIN_A];


y_CW_MIN_I= min(min(err_CW_I));
y_SL_MIN_I= min(min(err_SL_I));
y_LR_MIN_I= min(min(err_LR_I));
y_MIN_I = [y_CW_MIN_I y_SL_MIN_I  y_LR_MIN_I];

y_min= [y_MIN_P; y_MIN_A; y_MIN_I];
bar(x,y_min)
title('MIN Dimensional Error [km]')
legend('CW', 'SL' , 'LR')
set(gca,'xticklabel',name)

% MEAN ERROR
figure(3)
y_CW_MEAN_P = mean2(err_CW_P);
y_SL_MEAN_P = mean2(err_SL_P);
y_LR_MEAN_P = mean2(err_LR_P);
y_MEAN_P = [y_CW_MEAN_P y_SL_MEAN_P  y_LR_MEAN_P];

y_CW_MEAN_A = mean2(err_CW_A);
y_SL_MEAN_A = mean2(err_SL_A);
y_LR_MEAN_A = mean2(err_LR_A);
y_MEAN_A = [y_CW_MEAN_A y_SL_MEAN_A  y_LR_MEAN_A];


y_CW_MEAN_I = mean2(err_CW_I);
y_SL_MEAN_I = mean2(err_SL_I);
y_LR_MEAN_I = mean2(err_LR_I);
y_MEAN_I = [y_CW_MEAN_I y_SL_MEAN_I  y_LR_MEAN_I];

y_mean = [ y_MEAN_P; y_MEAN_A; y_MEAN_I];
bar(x,y_mean)
title('MEAN Dimensional Error [km]')
legend('CW', 'SL' , 'LR')
set(gca,'xticklabel',name)

% CC_plot = -d_theta:step:d_theta;
% TT_plot = TOF/3600;
% figure
% surf(CC_plot,TT_plot,err_CW_P','edgecolor', 'none');
% colorbar
% xlabel('\Delta \theta [°]')
% ylabel('TOF [h]')
% zlabel('C-W dimensional error [km]')


