function [] = plot_rdv_DRO(choice, choice_help, sol, Lam, Lam_t, r_f_syn, np, cr3bp, dro_T, dro_C, x0_T, x0_C, n_transfer, t_f, t_C)
%% Variables

p6 = zeros(1,np-2);
str = cell(1,np-2);

 
%% plots for choice = [1 0 0 0] (C-W)
 
if choice(1) && ~choice(2) && ~choice_help
    
   err_CW  = sol.CW.err1;
   fg_CW   = sol.CW.fg;
   
%  first guesses and Lambert solution synodic 
   figure(2)
   p3 = plot3(fg_CW(:,1), fg_CW(:,2), fg_CW(:,3),'blue');
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
          nom = ['First guess' 'Lambert opt' 'Start' str{:} 'Docking'];
          legend(pp,nom)
   else
      legend([p3 p4 p5 p7],'First guess', 'Lambert', 'Start', 'Docking')
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
    
   err_SL  = sol.SL.err1;
   fg_SL   = sol.SL.fg;
   
%  first guesses and Lambert solution synodic 
   figure(4)       
   p3 = plot3(fg_SL(:,1), fg_SL(:,2), fg_SL(:,3),'blue');
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
          nom = ['First guess' 'Lambert opt' 'Start' str{:} 'Docking'];
          legend(pp,nom)
   else
      legend([p3 p4 p5 p7],'First guess', 'Lambert', 'Start', 'Docking')
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
    
   err_LR  = sol.LR.err1;
   fg_LR   = sol.LR.fg;
   
   %  first guesses and Lambert solution synodic 
   figure(6)     
   p3 = plot3(fg_LR(:,1), fg_LR(:,2), fg_LR(:,3),'blue');
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
          nom = ['First guess' 'Lambert opt' 'Start' str{:} 'Docking'];
          legend(pp,nom)
   else
      legend([p3 p4 p5 p7],'First guess', 'Lambert', 'Start', 'Docking')
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
 
%% plots for choice = [1 1 1 0] (hybrid solution)

if choice(1) && choice(2)  && choice(3) && ~choice_help
    
   err_CW  = sol.CW.err1;
   err_SL  = sol.SL.err1;
   err_LR  = sol.LR.err1;
   fg_final = sol.final.fg;
   
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

%  first guesses and the best solution synodic 
   figure(6)     
   p3 = plot3(fg_final(:,1), fg_final(:,2), fg_final(:,3),'blue','LineWidth',3);
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
          nom = ['First guess' 'Best opt' 'Start' str{:} 'Docking'];
          legend(pp,nom)
   else
      legend([p3 p4 p5 p7],'First guess', 'Best', 'Start', 'Docking')
   end
 
   title('Luquette and the best Solution, Synodic Frame')
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
    
   err_cont = sol.cont.err(1);
   fg_cont  = sol.cont.fg;
   
%  first guesses and Lambert solution synodic 
   figure(8)
   p3 = plot3(fg_cont(:,1), fg_cont(:,2), fg_cont(:,3),'blue');
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

r = state_time(t_f(end),dro_T);
w = r(1);
e = r(2);
t = r(3);


% RDV in synodic frame
figure(10)
p1 = plot(dro_T.yv(:,1),dro_T.yv(:,2),'g');  %Target 
hold on
axis equal
grid on

p2 = plot(dro_C.yv(:,1),dro_C.yv(:,2),'k');  %Chaser 

% first guess
if choice(1) && sum(choice) == 1 % CW
    p3 = plot(fg_CW(:,1) + x0_T(1), fg_CW(:,2) + x0_T(2), 'blue'); 
elseif choice(2) && sum(choice) == 1 % SL
    p3 = plot(fg_SL(:,1) + x0_T(1), fg_SL(:,2) + x0_T(2), 'blue');
elseif choice(3) && sum(choice) == 1 % LR
    p3 = plot(fg_LR(:,1) + x0_T(1), fg_LR(:,2) + x0_T(2), 'blue');
elseif choice(4) && sum(choice) == 1 % continuation
    p3 = plot(fg_cont(:,1) + x0_T(1), fg_cont(:,2) + x0_T(2), 'blue');
else % [1 1 1 0]
    p3 = plot(fg_final(:,1) + x0_T(1), fg_final(:,2) + x0_T(2), 'blue');
end

p4 = plot(x0_T(1), x0_T(2), 'mo'); 
p5 = plot(x0_C(1), x0_C(2), 'co','LineWidth',1.5); 
p7 = plot(Lam(:,7) + Lam(:,1), Lam(:,8) + Lam(:,2), 'r','LineWidth',1.5);
p8 = plot(w, e, 'ko','LineWidth',1.5);

   
R_Moon = 1737/cr3bp.L; 
MOON = imread('moon.jpg','jpg');
props.FaceColor='texture';
props.EdgeColor='none';
props.FaceLighting='phong';
props.Cdata = MOON;
Center = [1 - cr3bp.mu; 0; 0];
[XX, YY, ZZ] = ellipsoid(-Center(1),Center(2),Center(3),R_Moon,R_Moon,R_Moon,30);
surface(-XX, -YY, -ZZ, props);

if np > 2
   for k = 2:np-1
       uu = state_time(t_f(k),dro_T);
       a = uu(1);
       b = uu(2);
       c = uu(3);
       p6(k-1) = plot3(r_f_syn(k,1) + a, r_f_syn(k,2) + b, r_f_syn(k,3) + c,'bo');
       str{k-1} = cellstr(sprintf('%s %d','HP', k-1));
   end
   pp = [p1 p2 p3 p4 p5 p6 p7 p8];
   nom = ['Target Orbit' 'Chaser Orbit' 'First guess' 'Target initial position' 'Chaser initial position' 'RDV trajectory' str{:} 'Docking'];
   legend(pp,nom)
else
   legend([p1 p2 p3 p4 p5 p7 p8],'Target Orbit', 'Chaser Orbit', 'First guess', 'Target initial position', 'Chaser initial position', 'RDV trajectory' ,'Docking')
end

 title('Rendezvous in the Synodic Frame')
 xlabel('x_{syn}')
 ylabel('y_{syn}')
 zlabel('z_{syn}')
 

%% RDV in LVLH frame
% Plot 2D for DRO

r = state_time(t_f(end),dro_T);
r_lvlh = syn2lvlh(t_f(end), r', r',  cr3bp.mu); % final position equal to target position: [0 0 0 0 0 0]'
w = r_lvlh(1);
e = r_lvlh(2);
t = r_lvlh(3);

 
 for i = 1:length(dro_C.yv)
    r = state_time(dro_C.tv(i), dro_T);
    yv_C_lvlh(1:6, i) = syn2lvlh(dro_C.tv(i), dro_C.yv(i,1:6)', r',  cr3bp.mu); % to convert the chaser state values in LVLH frame it's necessary to take the target state and the chaser state at the same time   
 end
 
 x0_C_lvlh = syn2lvlh(t_C, x0_C, x0_T,  cr3bp.mu); %Initial chaser position = initial point

for j = 1:length(Lam_t)
    Lam_lvlh(1:6, j) = rsyn2rlvlh(Lam_t(j), Lam(j,1:6)', Lam(j,7:12)',  cr3bp.mu);
end
 
 
 figure(12)
 p1 = plot(x0_C_lvlh(1),  x0_C_lvlh(3),'co','LineWidth',1.5);
 hold on
 axis equal
 grid on
 p2 = plot(yv_C_lvlh(1,:),yv_C_lvlh(3,:),'k');  %Chaser 
 p3 =  plot(Lam_lvlh(1,:), Lam_lvlh(3,:),'r','LineWidth',1.5); %RDV trajectory
 p7 = plot(w, t,'ko','LineWidth',1.5);

 if np > 2
 for k = 2:np-1
       uu = state_time(t_f(k),dro_T);
       uu_lvlh = syn2lvlh(t_f(k), uu', uu', cr3bp.mu);
       a = uu_lvlh(1);
       b = uu_lvlh(2);
       c = uu_lvlh(3);
       p4(k-1) = plot(r_f(k,1) + a, r_f(k,3) + c,'bo');
       str{k-1} = cellstr(sprintf('%s %d','HP', k-1));
 end
  pp = [p1 p2 p3 p4 p7];
  nom = ['Chaser initial position' 'Chaser Orbit' 'RDV trajectory' str{:} 'Docking' ];
  legend(pp,nom)
 else
   legend([p1 p2 p3 p7],'Chaser initial position', 'Chaser Orbit', 'RDV trajectory', 'Docking')
end
 
   
 title('Rendezvous in the LVLH Frame')
 xlabel('Downrange')
 ylabel('Altitude')
 
end