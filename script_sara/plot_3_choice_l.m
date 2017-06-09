function [] = plot_3_choice_l(choice, choice_help, sol, Lam1, Lam_t, r_f_syn, r_f , np, cr3bp, nro_T, nro_C, x0_T, x0_C, n_transfer, t_f,t_C, fg, err1)
%% Variables
Lam  = Lam1;

if choice(1) && ~choice(2) && ~choice_help
   err_CW  = sol.CW.err1;
   fg_CW   = sol.CW.fg;
end

if choice(2) && ~choice(1) && ~choice_help
   err_SL  = sol.SL.err1;
   fg_SL   = sol.SL.fg;
end

if choice(3) && ~choice(2) && ~choice_help
   err_LR  = sol.LR.err1;
   fg_LR   = sol.LR.fg;
end

if choice(1) && choice(2) && choice(3) && ~choice_help
   err_CW  = sol.CW.err1;
   err_SL  = sol.SL.err1;
   err_LR  = sol.LR.err1;
end

if choice(4) || choice_help
   err_cont = sol.cont.err(1);
   fg_cont  = sol.cont.fg;
end

p6 = zeros(1,np-2);
str = cell(1,np-2);

%% plots for choice = [1 1 1 0]

if choice(1) && choice(2)  && choice(3) && ~choice_help
%  Dimensional Error
   figure(1)
   semilogy(n_transfer,err_CW*cr3bp.L,'ro')
   hold on
   grid on
   semilogy(n_transfer,err_SL*cr3bp.L,'mo')
   semilogy(n_transfer,err_LR*cr3bp.L,'ko')
   
   set(gca,'XTick', 1:n_transfer);
   legend('CW_{syn}','SL_{syn}','L_{syn}')
   title('Init Dim Error for each Transfer, Syn Dynamics')
 
   xlabel('Transfer')
   ylabel('Error [km]')
 end
 
%% plots for choice = [1 0 0 0] (C-W)
 
if choice(1) && ~choice(2) && ~choice_help
%  first guesses and Lambert solution synodic 
   figure(2)
   p1 = plot3(fg_CW(:,1), fg_CW(:,2), fg_CW(:,3),'r');
   hold on
   grid on
   p4 = plot3(Lam(:,1), Lam(:,2), Lam(:,3),'g');
   p5 = plot3(r_f_syn(1,1), r_f_syn(1,2), r_f_syn(1,3),'co');
   p7 = plot3(0,0,0,'ko');
   
   if np > 2
      for k = 2:np-1
          p6(k-1) = plot3(r_f_syn(k,1), r_f_syn(k,2), r_f_syn(k,3),'bo');
          str{k-1} = cellstr(sprintf('%s %d','HP', k-1));
      end
          pp = [p1 p4 p5 p6 p7];
          nom = ['C-W' 'Lambert opt' 'Start' str{:} 'Docking'];
          legend(pp,nom)
   else
      legend([p1 p4 p5 p7],'C-W', 'Lambert', 'Start', 'Docking')
   end
 
 
   title('CW and Lambert Solution, Synodic Frame')
   xlabel('x_{syn}') 
   ylabel('y_{syn}')
   zlabel('z_{syn}')
 
 % Dimensional Error
   figure(3)
   semilogy(n_transfer, err_CW*cr3bp.L,'ro')
   hold on
   grid on

   set(gca,'XTick', 1:n_transfer);
   legend('CW_{syn}')
   title('Initial Dimensional Error for each Transfer, Synodic Dynamics')
 
   xlabel('Transfer')
   ylabel('Error [km]')
 
end
%% plots for choice = [0 1 0 0] (SL)


if choice(2)  && ~choice(1) && ~choice_help
%  first guesses and Lambert solution synodic 
   figure(4)       
   p2 = plot3(fg_SL(:,1), fg_SL(:,2), fg_SL(:,3),'m');
   hold on
   grid on
   p4 = plot3(Lam(:,1), Lam(:,2), Lam(:,3),'g');
   p5 = plot3(r_f_syn(1,1), r_f_syn(1,2), r_f_syn(1,3),'co');
   p7 = plot3(0,0,0,'ko');
   
   if np > 2
      for k = 2:np-1
          p6(k-1) = plot3(r_f_syn(k,1), r_f_syn(k,2), r_f_syn(k,3),'bo');
          str{k-1} = cellstr(sprintf('%s %d','HP', k-1));
      end
          pp = [p2 p4 p5 p6 p7];
          nom = ['SL' 'Lambert opt' 'Start' str{:} 'Docking'];
          legend(pp,nom)
   else
      legend([p2 p4 p5 p7],'SL', 'Lambert', 'Start', 'Docking')
   end
 
 
   title('SL and Lambert Solution, Synodic Frame')
   xlabel('x_{syn}') 
   ylabel('y_{syn}')
   zlabel('z_{syn}')
 
 % Dimensional Error
   figure(5)
   semilogy(n_transfer,err_SL*cr3bp.L,'mo')
   hold on
   grid on

   set(gca,'XTick', 1:n_transfer);
   legend('SL_{syn}')
   title('Initial Dimensional Error for each Transfer, Synodic Dynamics')
 
   xlabel('Transfer')
   ylabel('Error [km]')
end
 

%% plots for choice = [0 0 1 0] (Luquette)

if choice(3)  && ~choice(1) && ~choice_help
%  first guesses and Lambert solution synodic 
   figure(6)     
   p3 = plot3(fg_LR(:,1), fg_LR(:,2), fg_LR(:,3),'k');
   hold on
   grid on
   p4 = plot3(Lam(:,1), Lam(:,2), Lam(:,3),'g');
   p5 = plot3(r_f_syn(1,1), r_f_syn(1,2), r_f_syn(1,3),'co');
   p7 = plot3(0,0,0,'ko');
  
   if np > 2
      for k = 2:np-1
          p6(k-1) = plot3(r_f_syn(k,1), r_f_syn(k,2), r_f_syn(k,3),'bo');
          str{k-1} = cellstr(sprintf('%s %d','HP', k-1));
      end
          pp = [p3 p4 p5 p6 p7];
          nom = ['L' 'Lambert opt' 'Start' str{:} 'Docking'];
          legend(pp,nom)
   else
      legend([p3 p4 p5 p7],'L', 'Lambert', 'Start', 'Docking')
   end
 
   title('Luquette and Lambert Solution, Synodic Frame')
   xlabel('x_{syn}') 
   ylabel('y_{syn}')
   zlabel('z_{syn}')
 
%  Dimensional Error
   figure(7)
   semilogy(n_transfer,err_LR*cr3bp.L,'mo')
   hold on
   grid on

   set(gca,'XTick', 1:n_transfer);
   legend('L_{syn}')
   title('Initial Dimensional Error for each Transfer, Synodic Dynamics')
 
   xlabel('Transfer')
   ylabel('Error [km]')
 
end
 

%% plots for choice = [0 0 0 1] (Continuation)

if choice(4) || choice_help
%  first guesses and Lambert solution synodic 
   figure(8)
   p3 = plot3(fg_cont(:,1), fg_cont(:,2), fg_cont(:,3),'k');
   hold on
   grid on
   p4 = plot3(Lam(:,1), Lam(:,2), Lam(:,3),'g');
   p5 = plot3(r_f_syn(1,1), r_f_syn(1,2), r_f_syn(1,3),'co');
   p7 = plot3(0,0,0,'ko');
   
   if np > 2
      for k = 2:np-1
          p6(k-1) = plot3(r_f_syn(k,1), r_f_syn(k,2), r_f_syn(k,3),'bo');
          str{k-1} = cellstr(sprintf('%s %d','HP', k-1));
      end
          pp = [p3 p4 p5 p6 p7];
          nom = ['Continuation' 'Lambert opt' 'Start' str{:} 'Docking'];
          legend(pp,nom)
   else
      legend([p3 p4 p5 p7],'Continuation', 'Lambert', 'Start', 'Docking')
   end
 
 
   title('Continuation and Lambert Solution, Synodic Frame')
   xlabel('x_{syn}') 
   ylabel('y_{syn}')
   zlabel('z_{syn}')
 
%  Dimensional Error
   figure(9)
   semilogy(n_transfer, err_cont*cr3bp.L,'mo')
   hold on
   grid on
   
   set(gca,'XTick', 1:n_transfer);
   legend('Continuation_{syn}')
   title('Initial Dimensional Error for each Transfer, Synodic Dynamics')
 
   xlabel('Transfer')
   ylabel('Error [km]')
 
end

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
rad_sphere = 2/cr3bp.L;  % radius of the Approach Sphere to ensure trajectory safety policies [km]
p9 = sphere_plot(rad_sphere,w, e,t) % Approach Sphere centered on the Target
   
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
   for k = 2:np-1
       uu = state_time(t_f(k),nro_T);
       a = uu(1);
       b = uu(2);
       c = uu(3);
       p6(k-1) = plot3(r_f_syn(k,1) + a, r_f_syn(k,2) + b, r_f_syn(k,3) + c,'bo');
       str{k-1} = cellstr(sprintf('%s %d','HP', k-1));
   end
   pp = [p1 p2 p3 p4 p5 p6 p7 p9];
   nom = ['Target Orbit' 'Chaser Orbit' 'Target initial position' 'Chaser initial position' 'RDV trajectory' str{:} 'Docking' 'Approach Sphere'];
   legend(pp,nom)
else
   legend([p1 p2 p3 p4 p5 p7 p9],'Target Orbit', 'Chaser Orbit', 'Target initial position', 'Chaser initial position', 'RDV trajectory' ,'Docking', 'Approach Sphere' )
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
r_lvlh = syn2lvlh(t_f(end), r', r',  cr3bp.mu); % final position equal to target position: [0 0 0 0 0 0]'
w = r_lvlh(1);
e = r_lvlh(2);
t = r_lvlh(3);


 for j = 1:length(Lam)
   rT      = state_time(Lam_t(j), nro_T);
   rC      = state_time(Lam_t(j)- t_f(1) + t_C, nro_C);
   yv_C_lvlh(1:6, j) = syn2lvlh(Lam_t(j), rC', rT',  cr3bp.mu);
   
 end

 
 
 x0_C_lvlh = syn2lvlh(t_f(1), x0_C, x0_T,  cr3bp.mu); %Initial chaser position = initial point
 p1 = plot3(x0_C_lvlh(1), x0_C_lvlh(2),  x0_C_lvlh(3),'co','LineWidth',1.5);
 
 hold on
 axis equal
 grid on
 p2 = plot3(yv_C_lvlh(1,:),yv_C_lvlh(2,:),yv_C_lvlh(3,:),'k');  %Chaser Orbit   


for j = 1:length(Lam)
    Lam_lvlh(1:6, j) = rsyn2rlvlh(Lam_t(j), Lam(j,1:6)', Lam(j,7:12)',  cr3bp.mu);
end

for r = 1:size(fg,1)
    fg_lvlh(1:6, r) = rsyn2rlvlh(Lam_t(r), fg(r,1:6)', fg(r,7:12)',  cr3bp.mu);
end
 p3 = plot3(Lam_lvlh(1,:), Lam_lvlh(2,:), Lam_lvlh(3,:),'r','LineWidth',1.5); %RDV trajectory
 p7 = plot3(w, e, t,'ko','LineWidth',1.5);
 p8 = plot3(fg_lvlh(1,:), fg_lvlh(2,:), fg_lvlh(3, :),'b'); % First Guess
 p9 = sphere_plot(rad_sphere, 0, 0, 0); % Plot of the Approach Sphere in 3D centered on the Target
 
if np > 2
  for k = 2:np-1
       uu = state_time(t_f(k),nro_T);
       uu_lvlh = syn2lvlh(t_f(k), uu', uu', cr3bp.mu);
       a = uu_lvlh(1);
       b = uu_lvlh(2);
       c = uu_lvlh(3);
       p4(k-1) = plot3(r_f(k,1) + a, r_f(k,2) + b, r_f(k,3) + c,'bo');
       str{k-1} = cellstr(sprintf('%s %d','HP', k-1));
  end
   pp = [p1 p2 p3  p4 p7 p8 p9];
   nom = ['Chaser initial position' 'Chaser Orbit' 'RDV trajectory' str{:} 'Docking' 'First Guess' 'Approach Sphere'];
   legend(pp,nom)
  else
   legend([p1 p2 p3 p7 p8 p9],'Chaser initial position', 'Chaser Orbit', 'RDV trajectory', 'Docking', 'First Guess', 'Approach Sphere')
end

   
 title('Rendezvous in the LVLH Frame')
 xlabel('x_{LVLH}')
 ylabel('y_{LVLH}')
 zlabel('z_{LVLH}')
 
 % Plot 2D
 figure(12)
 p1=plot(x0_C_lvlh(1),  x0_C_lvlh(3),'co','LineWidth',1.5);
 hold on
 axis equal
 grid on
 p2 = plot(yv_C_lvlh(1,:),yv_C_lvlh(3,:),'k');  %Chaser 
 p3 =  plot(Lam_lvlh(1,:), Lam_lvlh(3,:),'r','LineWidth',1.5); %RDV trajectory
 p7 = plot(w, t,'ko','LineWidth',1.5);
 p8 = plot(fg_lvlh(1,:), fg_lvlh(3,:),'b'); % First Guess   
 p9 = circle(0,0,rad_sphere); % Plot of the Approach Sphere in 2D
 
 if np > 2
 for k = 2:np-1
       uu = state_time(t_f(k),nro_T);
       uu_lvlh = syn2lvlh(t_f(k), uu', uu', cr3bp.mu);
       a = uu_lvlh(1);
       b = uu_lvlh(2);
       c = uu_lvlh(3);
       p4(k-1) = plot(r_f(k,1) + a, r_f(k,3) + c,'bo');
%        p5(k-1) = plot(  a + err1,  c + err1,'ro');
       str{k-1} = cellstr(sprintf('%s %d','HP', k-1));
 end
  pp = [p1 p2 p3 p4 p7 p8 p9];
  nom = ['Chaser initial position' 'Chaser Orbit' 'RDV trajectory' str{:} 'Docking' 'First Guess' 'Approach Sphere'];
  legend(pp,nom)
 else
   legend([p1 p2 p3 p7 p8 p9],'Chaser initial position', 'Chaser Orbit', 'RDV trajectory', 'Docking', 'First Guess', 'Approach Sphere')
end
 
 title('Rendezvous in the LVLH Frame')
 xlabel('Downrange')
 ylabel('Altitude')
  

