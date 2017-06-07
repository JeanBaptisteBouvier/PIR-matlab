% Comparison between linear targeting algorithms
% 
% The first set of parameters allows to compare the errors when changing
% the radius of the DRO and the gap between the Chaser and Target orbits
% 
% The second set of parameters allows to compare the effect of the angular
% position on the orbit
% 

init;


% First set of parameters

gap = 2000; % gap between the chaser and target orbits (km)
size_inf = 16000; % radius of the smallest target orbit (km)
size_stp = 1000; % step between each target orbit (km)
size_sup = 300000; % radius of the biggest target orbit (km)



% Second set of parameters

% Ax = 240000; % radius of the target orbit (km)
% gap = 2000; % gap between the chaser and target orbits (km)
angl_inf = 0; % smallest angular position of the chaser (degree)
angl_stp = 2; % step between each angles (degree)
angl_sup = 346; % biggest angular position of the chaser (degree)
angl_diff = 10; % initial angle between target and chaser (degree)

Comp = struct('CW', [], 'SL', [], 'LR', []);

size_index = 1;
for Ax = size_inf:size_stp:size_sup
    sprintf(['Ax = ' num2str(Ax) ' km'])
    
    angl_index = 1;
    for theta = angl_inf:angl_stp:166


        %% Orbits of the Target and the Chaser in the synodic frame

        [dro_T, default, ~] = A_DRO(Ax, default, cst);
        theta_T = (theta + angl_diff)*pi/180;      % initial positions of the Target

        [dro_C, default, cr3bp] = A_DRO(Ax+gap, default, cst);
        theta_C = theta*pi/180;                   % initial positions of the Chaser


        %% Initialisation of the main structures
        hold_points_dim = []; % control point in the LVLH framework
        TOF_dim  = [3600*5];
        TOF_cont = [3600]; 
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

        output = Rendezvous_choice_dro(hold_points_dim, TOF_dim, init_cond, orbits, settings, continuation);

        Comp.CW(size_index, angl_index) = output.linear{1,1}.CW.err1;
        Comp.SL(size_index, angl_index) = output.linear{1,1}.SL.err1;
        Comp.LR(size_index, angl_index) = output.linear{1,1}.LR.err1;
        angl_index = angl_index+1;
    end

    for theta = 182:angl_stp:angl_sup


        %% Orbits of the Target and the Chaser in the synodic frame

        [dro_T, default, ~] = A_DRO(Ax, default, cst);
        theta_T = (theta + angl_diff)*pi/180;      % initial positions of the Target

        [dro_C, default, cr3bp] = A_DRO(Ax+gap, default, cst);
        theta_C = theta*pi/180;                   % initial positions of the Chaser


        %% Initialisation of the main structures
        hold_points_dim = []; % control point in the LVLH framework
        TOF_dim  = [3600*5];
        TOF_cont = [3600]; 
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

        output = Rendezvous_choice_dro(hold_points_dim, TOF_dim, init_cond, orbits, settings, continuation);

        Comp.CW(size_index, angl_index) = output.linear{1,1}.CW.err1;
        Comp.SL(size_index, angl_index) = output.linear{1,1}.SL.err1;
        Comp.LR(size_index, angl_index) = output.linear{1,1}.LR.err1;
        angl_index = angl_index+1;
    end
    size_index = size_index+1;
    
end

save Comp Comp

% figure;
% hold on
% grid on

Ax = size_inf:size_stp:size_sup ;
% plot(Ax, CW, 'red')
% % plot(Ax, SL, 'green')
% plot(Ax, LR, 'blue')
% % legend('Clohessy-Wiltshire', 'Straight-Line', 'Linearized Relative')
% legend('Clohessy-Wiltshire', 'Linearized Relative')
% title(['Error comparison between Targeting algorithms Dx = ' num2str(gap) ' km'])
% xlabel('Radius of the smaller DRO')


% plot(Theta, CW, 'red')
% plot(Theta, SL, 'green')
% plot(Theta, LR, 'blue')
% legend('Clohessy-Wiltshire', 'Straight-Line', 'Linearized Relative')
% title(['Error comparison between Targeting algorithms Dx = ' num2str(gap) ' km, Ax = ' num2str(Ax) ' km'])
% xlabel(['Initial angular position of the Chaser, ' num2str(angl_diff) ' degrees behind the Target'])
% 
% ylabel('error')

Theta = angl_inf:angl_stp:angl_sup ;
Theta_cut = horzcat(angl_inf:angl_stp:165, 182.5:angl_stp:angl_sup);
% Use of interpolation to find the points around 180°
l = length(Comp.CW(:,1));
L = length(Theta);
for k = 1:l
    CW(k,:) = interp1(Theta_cut, Comp.CW(k,:), Theta, 'spline');
    SL(k,:) = interp1(Theta_cut, Comp.SL(k,:), Theta, 'spline');
    LR(k,:) = interp1(Theta_cut, Comp.LR(k,:), Theta, 'spline');
end

figure;
hold on
grid on
[Theta,Size] = meshgrid(Theta, Ax);

surf(Size,Theta,CW, ones(l,L));
surf(Size,Theta,SL, zeros(l,L));
surf(Size,Theta,LR, 0.5*ones(l,L)); 

legend('Linearized Relative', 'Straight-Line', 'Clohessy-Wiltshire')
title('Error comparison between Targeting algorithms')
xlabel('Radius of the Target DRO, which is  2 000 km less than the radius of the Chaser DRO')
ylabel('Initial angular position of the Chaser, 10 degrees behind the Target')
zlabel('Error')

