function [] = plot_choice(Lam1, Lam_t, r_f_syn, r_f , np, cr3bp, nro_T, nro_C, x0_T, x0_C, t_f,t_C, primary ,disp_HP, color, trajectory)
%% Selection of LVLH reference primary

if(primary == 1)  % It will be choosen the ECI 
    syn2lvlh_c          = @syn2lvlh_e;
    rsyn2rlvlh_c        = @rsyn2rlvlh_e;
else              % It will be choosen the LCI 
    syn2lvlh_c          = @syn2lvlh;
    rsyn2rlvlh_c        = @rsyn2rlvlh;
end
%%
p6 = zeros(1,np-2);
str = cell(1,np-2);

 if disp_HP ==0
         l = 2:np-1;
 else
        l = 2:disp_HP;
        np = disp_HP;
 end
%% Variables
Lam  = Lam1;

%% plots plotted in any case

r = state_time(t_f(end),nro_T);
w = r(1);
e = r(2);
t = r(3);

% RDV in synodic frame
figure(10)
p1 = plot3(nro_T.yv(:,1),nro_T.yv(:,2),nro_T.yv(:,3),'g');  %Target 
hold on
axis equal
grid on
p2 = plot3(nro_C.yv(:,1),nro_C.yv(:,2),nro_C.yv(:,3),'k');  %Chaser 
p3 = plot3(x0_T(1), x0_T(2), x0_T(3),'mo');
p4 = plot3(x0_C(1), x0_C(2), x0_C(3),'co','LineWidth',1.5); 
p5 = plot3(Lam(:,7) + Lam(:,1), Lam(:,8) + Lam(:,2), Lam(:,9) + Lam(:,3),'r','LineWidth',1.5);
p7 = plot3(w, e, t,'ko','LineWidth',1.5);
rad_sphere = 2/cr3bp.L;  % radius of the Approach Sphere (2 km) to ensure trajectory safety policies [km]
rad_sphere_KOS = 200e-3/cr3bp.L; % radius of the Keep Out Sphere (200 m) to ensure trajectory safety policies [km]
p9 = sphere_plot(rad_sphere,w, e,t); % Approach Sphere centered on the Target
p10 = sphere_plot(rad_sphere_KOS, w, e, t); % Plot of the Keep Out Sphere in 3D centered on the Target

   
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
   for k = l
       uu = state_time(t_f(k),nro_T);
       a = uu(1);
       b = uu(2);
       c = uu(3);
       p6(k-1) = plot3(r_f_syn(k,1) + a, r_f_syn(k,2) + b, r_f_syn(k,3) + c,'bo');
       str{k-1} = cellstr(sprintf('%s %d','HP', k-1));
   end
   pp = [p1 p2 p3 p4 p5 p6 p7 p9 p10];
   nom = ['Target Orbit' 'Chaser Orbit' 'Target initial position' 'Chaser initial position' 'RDV trajectory' str{:} 'Docking' 'Approach Sphere' 'Keep Out Sphere'];
%    legend(pp,nom)
else
   legend([p1 p2 p3 p4 p5 p7 p9 p10],'Target Orbit', 'Chaser Orbit', 'Target initial position', 'Chaser initial position', 'RDV trajectory' ,'Docking', 'Approach Sphere', 'Keep Out Sphere' )
end

 title('Rendezvous in the Synodic Frame')
 xlabel('x_{syn}')
 ylabel('y_{syn}')
 zlabel('z_{syn}')
 
 %%
 
 % RDV in LVLH frame
 % Plot 3D
 
 figure(11)

r      = state_time(t_f(end),nro_T);
r_lvlh = syn2lvlh_c(t_f(end), r', r',  cr3bp.mu); % final position equal to target position: [0 0 0 0 0 0]'
w = r_lvlh(1);
e = r_lvlh(2);
t = r_lvlh(3);


 for j = 1:size(Lam,1)
   rT      = state_time(Lam_t(j), nro_T);
   rC      = state_time(Lam_t(j)- t_f(1) + t_C, nro_C);
   yv_C_lvlh(1:6, j) = syn2lvlh_c(Lam_t(j), rC', rT',  cr3bp.mu);
   
 end

 x0_C_lvlh = syn2lvlh_c(t_f(1), x0_C, x0_T,  cr3bp.mu); %Initial chaser position = initial point
 p1 = plot3(cr3bp.L*x0_C_lvlh(1), cr3bp.L*x0_C_lvlh(2),  cr3bp.L*x0_C_lvlh(3),'go','LineWidth',1.5);
 
 hold on
 axis equal
 grid on
 p2 = plot3(cr3bp.L*yv_C_lvlh(1,:),cr3bp.L*yv_C_lvlh(2,:),cr3bp.L*yv_C_lvlh(3,:),'k','LineWidth',0.5);  %Chaser Orbit   


for j = 1:size(Lam,1)
    Lam_lvlh(1:6, j) = rsyn2rlvlh_c(Lam_t(j), Lam(j,1:6)', Lam(j,7:12)',  cr3bp.mu);
end

 p3 = plot3(cr3bp.L*Lam_lvlh(1,:), cr3bp.L*Lam_lvlh(2,:), cr3bp.L*Lam_lvlh(3,:),color,'LineWidth',1.5); %RDV trajectory
 p7 = plot3(cr3bp.L*w, cr3bp.L*e, cr3bp.L*t,'ko','LineWidth',1.5);
 p9 = sphere_plot(cr3bp.L*rad_sphere, 0, 0, 0); % Plot of the Approach Sphere in 3D centered on the Target
 p10 = sphere_plot(cr3bp.L*rad_sphere_KOS, 0, 0, 0); % Plot of the Keep Out Sphere in 3D centered on the Target

 if np > 2
  for k = l
       uu = state_time(t_f(k),nro_T);
       uu_lvlh = syn2lvlh_c(t_f(k), uu', uu', cr3bp.mu);
       a = uu_lvlh(1);
       b = uu_lvlh(2);
       c = uu_lvlh(3);
       p4(k-1) = plot3(cr3bp.L*r_f(k,1) + cr3bp.L*a, cr3bp.L*r_f(k,2) + cr3bp.L*b, cr3bp.L*r_f(k,3) + cr3bp.L*c,'mo', 'LineWidth',1);
       str{k-1} = cellstr(sprintf('%s %d','HP', k-1));
  end
   pp = [p1 p2 p3 p4 p7 p9 p10];
   nom = ['Chaser initial position' 'Chaser Orbit' trajectory str{:} 'Docking' 'Approach Sphere' 'Keep Out Sphere'];
   legend(pp,nom)
else
   legend([p1 p2 p3 p7 p9 p10],'Chaser initial position', 'Chaser Orbit', trajectory, 'Docking', 'Approach Sphere', 'Keep Out Sphere')
end

 title('Rendezvous in the LVLH Frame')
 xlabel('x_{LVLH} [km]')
 ylabel('y_{LVLH} [km]')
 zlabel('z_{LVLH} [km]')
 
 % Plot 2D
 figure(12)
 p1=plot(cr3bp.L*x0_C_lvlh(1),  cr3bp.L*x0_C_lvlh(3),'go','LineWidth',1.5);
 hold on
 axis equal
 grid on
 p2 = plot(cr3bp.L*yv_C_lvlh(1,:),cr3bp.L*yv_C_lvlh(3,:),'k','LineWidth',1.5);  %Chaser 
 p3 =  plot(cr3bp.L*Lam_lvlh(1,:), cr3bp.L*Lam_lvlh(3,:),color,'LineWidth',1.5); %RDV trajectory
 p7 = plot(cr3bp.L*w, cr3bp.L*t,'ko','LineWidth',1.5);
 p9 = circle(0,0,cr3bp.L*rad_sphere,'g'); % Plot of the Approach Sphere in 2D
 p10 = circle(0,0,cr3bp.L*rad_sphere_KOS, 'm'); % Plot of the Keep out Sphere in 2D
 
 if np > 2
 for k = l
       uu = state_time(t_f(k),nro_T);
       uu_lvlh = syn2lvlh_c(t_f(k), uu', uu', cr3bp.mu);
       a = uu_lvlh(1);
       b = uu_lvlh(2);
       c = uu_lvlh(3);
       p4(k-1) = plot(cr3bp.L*r_f(k,1) + cr3bp.L*a, cr3bp.L*r_f(k,3) + cr3bp.L*c,'mo','LineWidth',1);
       str{k-1} = cellstr(sprintf('%s %d','HP', k-1));
 end
  pp = [p1 p2 p3 p4 p7 p9 p10];
  nom = ['Chaser initial position' 'Chaser Orbit' trajectory str{:} 'Docking' 'Approach Sphere' 'Keep out Sphere'];
  legend(pp,nom)
 else
   legend([p1 p2 p3 p7 p9 p10],'Chaser initial position', 'Chaser Orbit', trajectory, 'Docking', 'Approach Sphere', 'Keep out Sphere')
end
 
 title('Rendezvous in the LVLH Frame')
 xlabel('Downrange [km]')
 ylabel('Altitude [km]')
  

