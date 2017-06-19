function [first_guess, Lam_fin, Lam_t, it, errv] = Lambert_cycle(x0_rel, r_next, x0_T, phi0, TOF, toll, it_max, cr3bp, options, tt)

% Variables Intialization
%Lam(it_max+1) = struct();   
%errv          = zeros(it_max+1,1);
err           = 1000;
it            = -1;  
i             = 1;
 
while err > toll && it < it_max
    % MATLAB slow computation, working on Windows and Linux 
%     [Lam_t, Lam] = ode113(@(t,y)crtbp(t,y,cr3bp.mu),[0 TOF],[x0_rel, x0_T, phi0']', options); 

    % MEX fast computation, on LINUX only
    [~, ~, Lam_t, Lam] = ode78_cr3bp_rel([0 TOF], [x0_rel, x0_T, phi0']', cr3bp.mu);
    Lam = Lam_shape(Lam); % change lines in columns
    
    Phi_Lam = reshape(Lam(end,13:48),6,6);
    Phi_rv = Phi_Lam(1:3,4:6); 
    
    r_obt = Lam(end,1:3);
    r_desired = r_next; 
    
    err = norm(r_obt - r_desired);
    
    dr0 = (r_obt - r_desired)';    
    dV_Lam = (Phi_rv) \ (dr0);
    x0_rel(1,4:6) = x0_rel(1,4:6) - dV_Lam';
    
    errv(i) = err;
    Lam_state{i} = Lam;
    
    i=i+1;
    it=it+1;
    
end

first_guess = Lam_state{1};
Lam_fin = Lam_state{end};





%  figure(7)
%  hold on
%  plot(Lam(i-1).yv(:,8)+Lam(i-1).yv(:,2), Lam(i-1).yv(:,9)+Lam(i-1).yv(:,3),'b','LineWidth',1.5);
