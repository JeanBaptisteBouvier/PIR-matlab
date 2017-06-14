function [] = plot_rdv_DRO(choice, sol, r_f_syn, r_f, np, cr3bp, dro_T, dro_C, x0_T, x0_C, t_f, t_C)
%% Variables




p6 = zeros(1,np-2);
str = cell(1,np-2);

%% First guesses and the solution synodic

figure(1)
hold on
grid on
Lam = sol{1}.Lam;
fg = sol{1}.fg;
p3 = plot3(fg(:,1), fg(:,2), fg(:,3),'b');
p4 = plot3(Lam(:,1), Lam(:,2), Lam(:,3),'g');
for i=2:np-1
   Lam = sol{i}.Lam;
   fg = sol{i}.fg;
   
   p3 = plot3(fg(:,1), fg(:,2), fg(:,3),'b');
   p4 = plot3(Lam(:,1), Lam(:,2), Lam(:,3),'g');
   p6(i-1) = plot3(r_f_syn(i,1), r_f_syn(i,2), r_f_syn(i,3),'bo');
   str{i-1} = cellstr(sprintf('%s %d','HP', i-1));
end
p5 = plot3(r_f_syn(1,1), r_f_syn(1,2), r_f_syn(1,3),'co');
p7 = plot3(0,0,0,'ko');
pp = [p3 p4 p5 p6 p7];
nom = ['First guess' 'Lambert opt' 'Start' str{:} 'Docking'];
legend(pp,nom)

if choice(1) && sum(choice) == 1
   title('CW and Lambert Solution, Synodic Frame')
elseif choice(2) && sum(choice) == 1
    title('SL and Lambert Solution, Synodic Frame')
elseif choice(3) && sum(choice) == 1
    title('Luquette and Lambert Solution, Synodic Frame')   
elseif sum(choice) == 3
    title('Best Solution and Lambert Solution, Synodic Frame')
end
xlabel('x_{syn}') 
ylabel('y_{syn}')
zlabel('z_{syn}')
hold off

%% Dimensional Error

figure(2)
hold on
grid on
if sum(choice) == 1
    for i=1:np-1
       err = sol{i}.err1;
       semilogy(i, err*cr3bp.L,'ro')
    end
    set(gca,'XTick', 1:np-1);
    
    if choice(1)
       legend('CW_{syn}')
    elseif choice(2)
       legend('SL_{syn}')
    elseif choice(3)
        legend('L_{syn}')
    end

elseif sum(choice) == 3
   for i=1:np-1
        err_CW  = sol{i}.err_CW;
        err_SL  = sol{i}.err_SL;
        err_LR  = sol{i}.err_LR;
        semilogy(i,err_CW*cr3bp.L,'ro')
        semilogy(i,err_SL*cr3bp.L,'mo')
        semilogy(i,err_LR*cr3bp.L,'ko')
    end
   set(gca,'XTick', 1:np-1);
   legend('CW_{syn}','SL_{syn}','L_{syn}')
end

title('Init Dim Error for each Transfer, Syn Dynamics')
xlabel('Number of the hold point')
ylabel('Error [km]')
 


%% plots plotted in any case

r = state_time(t_f(end),dro_T);
w = r(1);
e = r(2);


% RDV in synodic frame
figure(3)
hold on
axis equal
grid on

R_Moon = 1737/cr3bp.L; 
MOON = imread('moon.jpg','jpg');
props.FaceColor='texture';
props.EdgeColor='none';
props.FaceLighting='phong';
props.Cdata = MOON;
Center = [1 - cr3bp.mu; 0; 0];
[XX, YY, ZZ] = ellipsoid(-Center(1),Center(2),Center(3),R_Moon,R_Moon,R_Moon,30);
surface(-XX, -YY, -ZZ, props);

p1 = plot(dro_T.yv(:,1),dro_T.yv(:,2),'g');              % orbit of the Target
p2 = plot(dro_C.yv(:,1),dro_C.yv(:,2),'k');              % orbit of the Chaser 
p3 = plot(fg(:,1) + x0_T(1), fg(:,2) + x0_T(2), 'b','LineWidth',2.5); % first guess
p4 = plot(x0_T(1), x0_T(2), 'mo');                       % Target initial position
p5 = plot(x0_C(1), x0_C(2), 'co');                       % Chaser initial position
Lam = sol{1}.Lam;
p6 = plot(Lam(:,7) + Lam(:,1), Lam(:,8) + Lam(:,2), 'r','LineWidth',1.5);  % rdv trajectory

for i = 2:np-1 %hold points except start and arrival
   uu = state_time(t_f(i),dro_T);
   a = uu(1);
   b = uu(2);
   Lam = sol{i}.Lam;
   p6 = plot(Lam(:,7) + Lam(:,1), Lam(:,8) + Lam(:,2), 'r','LineWidth',1.5);  % rdv trajectory
   p7(i-1) = plot(r_f_syn(i,1) + a, r_f_syn(i,2) + b,'bo');                    % hold points
   str{i-1} = cellstr(sprintf('%s %d','HP', i-1));
end
p8 = plot(w, e, 'ko','LineWidth',1.5);
pp = [p1 p2 p3 p4 p5 p6 p7 p8];
nom = ['Target Orbit' 'Chaser Orbit' 'First guess' 'Target initial position' 'Chaser initial position' 'RDV trajectory' str{:} 'Docking'];
legend(pp,nom)


title('Rendezvous in the Synodic Frame')
xlabel('x_{syn}')
ylabel('y_{syn}')
zlabel('z_{syn}')
 

%% RDV in LVLH frame
% Plot 2D for DRO
figure(4)
hold on
axis equal
grid on

r = state_time(t_f(end),dro_T);
r_lvlh = syn2lvlh(t_f(end), r', r',  cr3bp.mu); % final position equal to target position: [0 0 0 0 0 0]'
w = r_lvlh(1);
t = r_lvlh(3);
x0_C_lvlh = syn2lvlh(t_f(1), x0_C, x0_T,  cr3bp.mu); %Initial chaser position = initial point

p1 = plot(x0_C_lvlh(1),  x0_C_lvlh(3),'co','LineWidth',1.5);

p7 = plot(w, t,'ko','LineWidth',1.5);

Lam = sol{1}.Lam;
Lam_t = sol{1}.Lam_t;
fg = sol{1}.fg;

for j = 1:length(Lam_t)
    rT      = state_time(Lam_t(j), dro_T);
    rC      = state_time(Lam_t(j)- t_f(1) + t_C, dro_C);
    yv_C_lvlh(1:6, j) = syn2lvlh(Lam_t(j), rC', rT',  cr3bp.mu);
    Lam_lvlh(1:6, j) = rsyn2rlvlh(Lam_t(j), Lam(j,1:6)', Lam(j,7:12)',  cr3bp.mu);
end
% Due to iterations of ode113, fg and Lambert solution have not necessarily
% the same size
for r = 1:min(length(Lam_t), length(fg(:,1)))
    fg_lvlh(1:6, r) = rsyn2rlvlh(Lam_t(r), fg(r,1:6)', fg(r,7:12)',  cr3bp.mu);
end
p2 = plot(yv_C_lvlh(1,:),yv_C_lvlh(3,:),'k');     % Chaser
p3 = plot(Lam_lvlh(1,:), Lam_lvlh(3,:),'r','LineWidth',1.5); % RDV trajectory
p8 = plot(fg_lvlh(1,:), fg_lvlh(3, :),'b');       % First Guess

for i=2:np-1
    clear Lam_lvlh
    clear fg_lvlh
    Lam = sol{i}.Lam;
    Lam_t = sol{i}.Lam_t;
    fg = sol{i}.fg;
    for j = 1:length(Lam_t)
        Lam_lvlh(1:6, j) = rsyn2rlvlh(Lam_t(j), Lam(j,1:6)', Lam(j,7:12)',  cr3bp.mu);
    end
    for r = 1:min(length(Lam_t), length(fg(:,1)))
        fg_lvlh(1:6, r) = rsyn2rlvlh(Lam_t(r), fg(r,1:6)', fg(r,7:12)',  cr3bp.mu);
    end
    p3 = plot(Lam_lvlh(1,:), Lam_lvlh(3,:),'r','LineWidth',1.5); %RDV trajectory
    p8 = plot(fg_lvlh(1,:), fg_lvlh(3, :),'b'); % First Guess

    uu = state_time(t_f(i),dro_T);
    uu_lvlh = syn2lvlh(t_f(i), uu', uu', cr3bp.mu);
    a = uu_lvlh(1);
    c = uu_lvlh(3);
%     p6(i-1) = plot(r_f(i,1) + a, r_f(i,3) + c,'bo');
    p6(i-1) = plot(r_f(i,1), r_f(i,3), 'bo');
    str{i-1} = cellstr(sprintf('%s %d','HP', i-1));
   
end
pp = [p1 p2 p3 p6 p7 p8];
nom = ['Chaser initial position' 'Chaser Orbit' 'RDV trajectory' str{:} 'Docking' 'First Guess'];
legend(pp,nom)
   
title('Rendezvous in the LVLH Frame')
xlabel('Downrange')
ylabel('Altitude')
 
end