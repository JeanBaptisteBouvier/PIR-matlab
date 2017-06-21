function [] = plot_choice_fg(choice, sol, r_f_syn, r_f, np, cr3bp, nro_T, nro_C, x0_T, x0_C, t_f, t_C,primary)
%% Selection of LVLH reference primary

if(primary == 1)  % It will be choosen the ECI 
    syn2lvlh_c          = @syn2lvlh_e;
    rsyn2rlvlh_c        = @rsyn2rlvlh_e;
else              % It will be choosen the LCI 
    syn2lvlh_c          = @syn2lvlh;
    rsyn2rlvlh_c        = @rsyn2rlvlh;
end

%% Variables
p6 = zeros(1,np-2);
str = cell(1,np-2);
k = 1:np-1;
l = 2:np-1;

%% First guesses and the solution synodic

figure(1)
hold on
grid on
Lam = sol{1}.Lam;
fg = sol{1}.fg;
p3 = plot3(cr3bp.L*fg(:,1), cr3bp.L*fg(:,2), cr3bp.L*fg(:,3),'b');
p4 = plot3(cr3bp.L*Lam(:,1), cr3bp.L*Lam(:,2),cr3bp.L*Lam(:,3),'g');
p5 = plot3(cr3bp.L*r_f_syn(1,1), cr3bp.L*r_f_syn(1,2), cr3bp.L*r_f_syn(1,3),'co');
p7 = plot3(0,0,0,'ko');
if np > 2
    for i=l
       Lam = sol{i}.Lam;
       fg = sol{i}.fg;

       p3 = plot3(cr3bp.L*fg(:,1), cr3bp.L*fg(:,2), cr3bp.L*fg(:,3),'b');
       p4 = plot3(cr3bp.L*Lam(:,1), cr3bp.L*Lam(:,2), cr3bp.L*Lam(:,3),'g');
       p6(i-1) = plot3(cr3bp.L*r_f_syn(i,1), cr3bp.L*r_f_syn(i,2), cr3bp.L*r_f_syn(i,3),'bo');
       str{i-1} = cellstr(sprintf('%s %d','HP', i-1));
    end
    pp = [p3 p4 p5 p6 p7];
    nom = ['First guess' 'Lambert opt' 'Start' str{:} 'Docking'];
     legend(pp,nom)
else
   legend([p3 p4 p5 p7],'First, guess', 'Lambert opt' ,'Start', 'Docking')
end

if choice(1) && sum(choice) == 1
   title('CW and Lambert Solution, Synodic Frame')
elseif choice(2) && sum(choice) == 1
    title('SL and Lambert Solution, Synodic Frame')
elseif choice(3) && sum(choice) == 1
    title('LR and Lambert Solution, Synodic Frame')   
elseif sum(choice) == 3
    title('Best Solution and Lambert Solution, Synodic Frame')
end
xlabel('x_{syn [km]}') 
ylabel('y_{syn [km]}')
zlabel('z_{syn [km]}')
hold off

%% Dimensional Error

figure(2)
hold on
grid on
if sum(choice) == 1
    for i=k
       err = sol{i}.err1;
       semilogy(i, err*cr3bp.L,'ro')
    end
    set(gca,'XTick', k);
    
    if choice(1)
       legend('CW_{syn}')
    elseif choice(2)
       legend('SL_{syn}')
    elseif choice(3)
        legend('L_{syn}')
    end

elseif sum(choice) == 3
   for i=k
        err_CW  = sol{i}.err_CW;
        err_SL  = sol{i}.err_SL;
        err_LR  = sol{i}.err_LR;
        semilogy(i,err_CW*cr3bp.L,'ro')
        semilogy(i,err_SL*cr3bp.L,'mo')
        semilogy(i,err_LR*cr3bp.L,'ko')
    end
   set(gca,'XTick', k);
   legend('CW_{syn}','SL_{syn}','L_{syn}')
end

title('Init Dim Error for each Transfer, Syn Dynamics')
xlabel('Number of the hold point')
ylabel('Error [km]')
 


%% plots plotted in any case

r = state_time(t_f(end),nro_T);
w = r(1);
e = r(2);
t = r(3);


% RDV in synodic frame
figure(3)

p1 = plot3(cr3bp.L*nro_T.yv(:,1),cr3bp.L*nro_T.yv(:,2),cr3bp.L*nro_T.yv(:,3),'g');              % orbit of the Target
hold on
axis equal
grid on

p2 = plot3(cr3bp.L*nro_C.yv(:,1),cr3bp.L*nro_C.yv(:,2),cr3bp.L*nro_C.yv(:,3),'k');              % orbit of the Chaser 
p3 = plot3(cr3bp.L*fg(:,1) + cr3bp.L*x0_T(1), cr3bp.L*fg(:,2) + cr3bp.L*x0_T(2),cr3bp.L*fg(:,3) + cr3bp.L*x0_T(3), 'b','LineWidth',2.5); % first guess
p4 = plot3(cr3bp.L*x0_T(1), cr3bp.L*x0_T(2), cr3bp.L*x0_T(3), 'mo');                       % Target initial position
p5 = plot3(cr3bp.L*x0_C(1), cr3bp.L*x0_C(2), cr3bp.L*x0_C(3), 'co');                       % Chaser initial position
Lam = sol{1}.Lam;
p6 = plot3(cr3bp.L*Lam(:,7) + cr3bp.L*Lam(:,1), cr3bp.L*Lam(:,8) + cr3bp.L*Lam(:,2), cr3bp.L*Lam(:,9) + cr3bp.L*Lam(:,3), 'r','LineWidth',1.5);  % rdv trajectory
rad_sphere = 2/cr3bp.L*cr3bp.L;  % radius of the Approach Sphere to ensure trajectory safety policies [km]
p9 = sphere_plot(rad_sphere,cr3bp.L*w, cr3bp.L*e,cr3bp.L*t); % Approach Sphere centered on the Target
p8 = plot3(cr3bp.L*w, cr3bp.L*e, cr3bp.L*t, 'ko','LineWidth',1.5);

R_Moon = 1737/cr3bp.L; 
MOON = imread('moon.jpg','jpg');
props.FaceColor='texture';
props.EdgeColor='none';
props.FaceLighting='phong';
props.Cdata = MOON;
Center = [1 - cr3bp.mu; 0; 0];
[XX, YY, ZZ] = ellipsoid(-Center(1),Center(2),Center(3),R_Moon,R_Moon,R_Moon,30);
surface(-XX, -YY, -ZZ,props);

if np > 2
    for i = l %hold points except start and arrival
       uu = state_time(t_f(i),nro_T);
       a = uu(1);
       b = uu(2);
       c = uu(3);
       Lam = sol{i}.Lam;
       p6 = plot3(cr3bp.L*Lam(:,7) + cr3bp.L*Lam(:,1), cr3bp.L*Lam(:,8) + cr3bp.L*Lam(:,2), cr3bp.L*Lam(:,9) + cr3bp.L*Lam(:,3), 'r','LineWidth',1.5);  % rdv trajectory
       p7(i-1) = plot3(cr3bp.L*r_f_syn(i,1) + cr3bp.L*a, cr3bp.L*r_f_syn(i,2) + cr3bp.L*b, cr3bp.L*r_f_syn(i,3) + cr3bp.L*c, 'bo') ;                   % hold points
       str{i-1} = cellstr(sprintf('%s %d','HP', i-1));
    end
    pp = [p1 p2 p3 p4 p5 p6 p7 p8 p9];
    nom = ['Target Orbit' 'Chaser Orbit' 'First guess' 'Target initial position' 'Chaser initial position' 'RDV trajectory' str{:} 'Docking' 'Approach Sphere'];
    legend(pp,nom)
else
   legend([p1 p2 p3 p4 p5 p6 p8 p9],'Target Orbit', 'Chaser Orbit' ,'First guess', 'Target initial position', 'Chaser initial position' ,'RDV trajectory', 'Docking' ,'Approach Sphere')
end

title('Rendezvous in the Synodic Frame')
xlabel('x_{syn [km]}')
ylabel('y_{syn [km]}')
zlabel('z_{syn [km]}')
 
%% RDV in LVLH frame
% Plot 3D for NRO 
figure(4)

r = state_time(t_f(end),nro_T);
r_lvlh = syn2lvlh_c(t_f(end), r', r',  cr3bp.mu); % final position equal to target position: [0 0 0 0 0 0]'
w = r_lvlh(1);
e = r_lvlh(2);
t = r_lvlh(3);
x0_C_lvlh = syn2lvlh_c(t_f(1), x0_C, x0_T,  cr3bp.mu); %Initial chaser position = initial point
p1 = plot3(cr3bp.L*x0_C_lvlh(1), cr3bp.L*x0_C_lvlh(2),  cr3bp.L*x0_C_lvlh(3),'co','LineWidth',1.5);
hold on
axis equal
grid on


Lam = sol{1}.Lam;
Lam_t = sol{1}.Lam_t;
fg = sol{1}.fg;

for j = 1:size(Lam,1)
    rT      = state_time(Lam_t(j), nro_T);
    rC      = state_time(Lam_t(j)- t_f(1) + t_C, nro_C);
    yv_C_lvlh(1:6, j) = syn2lvlh_c(Lam_t(j), rC', rT',  cr3bp.mu);
    Lam_lvlh(1:6, j) = rsyn2rlvlh_c(Lam_t(j), Lam(j,1:6)', Lam(j,7:12)',  cr3bp.mu);
end
p3 = plot3(cr3bp.L*Lam_lvlh(1,:), cr3bp.L*Lam_lvlh(2,:), cr3bp.L*Lam_lvlh(3,:),'r','LineWidth',1.5); % RDV trajectory
      
p7 = plot3(cr3bp.L*w, cr3bp.L*e, cr3bp.L*t,'ko','LineWidth',1.5);

% Due to iterations of ode113, fg and Lambert solution have not necessarily
% the same size
for r = 1:min(length(Lam_t), length(fg(:,1)))
    fg_lvlh(1:6, r) = rsyn2rlvlh(Lam_t(r), fg(r,1:6)', fg(r,7:12)',  cr3bp.mu);
end
p2 = plot3(cr3bp.L*yv_C_lvlh(1,:), cr3bp.L*yv_C_lvlh(2,:),cr3bp.L*yv_C_lvlh(3,:),'k','LineWidth',1.5);     % Chaser
p9 = sphere_plot(rad_sphere, 0, 0, 0); % Plot of the Approach Sphere in 3D centered on the Target
p8 = plot3(cr3bp.L*fg_lvlh(1,:),cr3bp.L*fg_lvlh(2, :), cr3bp.L*fg_lvlh(3, :),'b','LineWidth',1);       % First Guess


for i=l
    clear Lam_lvlh
    clear fg_lvlh
    Lam = sol{i}.Lam;
    Lam_t = sol{i}.Lam_t;
    fg = sol{i}.fg;
    for j = 1:length(Lam_t)
        Lam_lvlh(1:6, j) = rsyn2rlvlh_c(Lam_t(j), Lam(j,1:6)', Lam(j,7:12)',  cr3bp.mu);
    end
    for r = 1:min(length(Lam_t), length(fg(:,1)))
        fg_lvlh(1:6, r) = rsyn2rlvlh(Lam_t(r), fg(r,1:6)', fg(r,7:12)',  cr3bp.mu);
    end
    p8 = plot3(cr3bp.L*fg_lvlh(1,:),cr3bp.L*fg_lvlh(2, :), cr3bp.L*fg_lvlh(3, :),'b','LineWidth',1); % First Guess
    p3 = plot3(cr3bp.L*Lam_lvlh(1,:),cr3bp.L*Lam_lvlh(2,:), cr3bp.L*Lam_lvlh(3,:),'r','LineWidth',1.5); %RDV trajectory
end


if np > 2
    for i=l
    uu = state_time(t_f(i),nro_T);
    uu_lvlh = syn2lvlh_c(t_f(i), uu', uu', cr3bp.mu);
    a = uu_lvlh(1);
    b = uu_lvlh(2);
    c = uu_lvlh(3);
    p6(i-1) = plot3(cr3bp.L*r_f(i,1)+cr3bp.L*a , cr3bp.L*r_f(i,2)+cr3bp.L*b , cr3bp.L*r_f(i,3)+cr3bp.L*c , 'bo','LineWidth',1);
    str{i-1} = cellstr(sprintf('%s %d','HP', i-1));
    end 
    pp = [p1 p2 p3 p6 p7 p8 p9];
    nom = ['Chaser initial position' 'Chaser Orbit' 'RDV trajectory' str{:} 'Docking' 'First Guess' 'Approach Sphere'];     
    legend(pp,nom)
else
    legend([p1 p2 p3 p7 p8 p9],'Chaser initial position', 'Chaser Orbit', 'RDV trajectory', 'Docking','First Guess','Approach Sphere')
end
    

title('Rendezvous in the LVLH Frame')
xlabel('Downrange [km]')
ylabel('Altitude [km]')



% Plot 2D for NRO
figure(5)
hold on
axis equal
grid on
p1 = plot(cr3bp.L*x0_C_lvlh(1),  cr3bp.L*x0_C_lvlh(3),'co','LineWidth',1.5);

Lam = sol{1}.Lam;
Lam_t = sol{1}.Lam_t;
fg = sol{1}.fg;

for j = 1:size(Lam,1)
    rT      = state_time(Lam_t(j), nro_T);
    rC      = state_time(Lam_t(j)- t_f(1) + t_C, nro_C);
    yv_C_lvlh(1:6, j) = syn2lvlh_c(Lam_t(j), rC', rT',  cr3bp.mu);
    Lam_lvlh(1:6, j) = rsyn2rlvlh_c(Lam_t(j), Lam(j,1:6)', Lam(j,7:12)',  cr3bp.mu);
end
p3 = plot(cr3bp.L*Lam_lvlh(1,:), cr3bp.L*Lam_lvlh(3,:),'r','LineWidth',1.5); % RDV trajectory


p7 = plot(cr3bp.L*w, cr3bp.L*t,'ko','LineWidth',1.5);
% Due to iterations of ode113, fg and Lambert solution have not necessarily
% the same size
for r = 1:min(length(Lam_t), length(fg(:,1)))
    fg_lvlh(1:6, r) = rsyn2rlvlh(Lam_t(r), fg(r,1:6)', fg(r,7:12)',  cr3bp.mu);
end
p2 = plot(cr3bp.L*yv_C_lvlh(1,:),cr3bp.L*yv_C_lvlh(3,:),'k','LineWidth',1.5);     % Chaser
p9 = circle(0,0,rad_sphere,'r'); % Plot of the Approach Sphere in 2D
p8 = plot(cr3bp.L*fg_lvlh(1,:), cr3bp.L*fg_lvlh(3, :),'b','LineWidth',1);       % First Guess
for i=l
    clear Lam_lvlh
    clear fg_lvlh
        Lam = sol{i}.Lam;
        Lam_t = sol{i}.Lam_t;
        fg = sol{i}.fg;
        for j = 1:length(Lam_t)
            Lam_lvlh(1:6, j) = rsyn2rlvlh_c(Lam_t(j), Lam(j,1:6)', Lam(j,7:12)',  cr3bp.mu);
        end
        for r = 1:min(length(Lam_t), length(fg(:,1)))
            fg_lvlh(1:6, r) = rsyn2rlvlh(Lam_t(r), fg(r,1:6)', fg(r,7:12)',  cr3bp.mu);
        end
        p3 = plot(cr3bp.L*Lam_lvlh(1,:),cr3bp.L* Lam_lvlh(3,:),'r','LineWidth',1.5); %RDV trajectory
        p8 = plot(cr3bp.L*fg_lvlh(1,:), cr3bp.L*fg_lvlh(3, :),'b','LineWidth',1 );% First Guess
end

if np > 2
   for i=l
    uu = state_time(t_f(i),nro_T);
    uu_lvlh = syn2lvlh_c(t_f(i), uu', uu', cr3bp.mu);
    a = uu_lvlh(1);
    c = uu_lvlh(3);
    p6(i-1) = plot(cr3bp.L*r_f(i,1), cr3bp.L*r_f(i,3), 'bo','LineWidth',1);
    str{i-1} = cellstr(sprintf('%s %d','HP', i-1));
   end
    pp = [p1 p2 p3 p6 p7 p8 p9];
    nom = ['Chaser initial position' 'Chaser Orbit' 'RDV trajectory' str{:} 'Docking' 'First Guess' 'Approach Sphere'];
    legend(pp,nom)
else
    legend([p1 p2 p3 p7 p8 p9],'Chaser initial position', 'Chaser Orbit' ,'RDV trajectory', 'Docking','First Guess', 'Approach Sphere')
end
   
title('Rendezvous in the LVLH Frame')
xlabel('Downrange [km]')
ylabel('Altitude [km]')
 
end