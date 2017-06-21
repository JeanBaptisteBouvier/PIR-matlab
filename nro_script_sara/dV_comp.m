function output = dV_comp(Lam, x0_rel, np, cr3bp, disp_HP, MC_runs, sol_HP)

% check the presence of the dispersion
if disp_HP ==0
        k = 1:np-1;
        j = 2:np-1;
else
        k = 1:disp_HP;
        j = 2:disp_HP;
end

% initial impulse from NRO 
dV1 = zeros(MC_runs+1, 1);
if  disp_HP == 1 
    dV1(1,1) = norm(Lam{1}(1,4:6) - x0_rel(4:6)');
    for l = 1: MC_runs
     Lam_HP{l} = sol_HP{1+l}.Lam_HP;
     dV1(l+1,1) = norm(Lam_HP{l}(1,4:6) - x0_rel(4:6)');
    end
  else
   dV1(1:MC_runs+1,1) = norm(Lam{1}(1,4:6) - x0_rel(4:6)');
end

% velocity variation for each hold point
dV2 = zeros(MC_runs+1, length(j));
for i = j
  if i == disp_HP
    dV2(1,i) = norm(Lam{i}(1,4:6) - Lam{i-1}(end ,4:6));
    for l = 1: MC_runs
     Lam_HP{l} = sol_HP{i+l}.Lam_HP;
     dV2(l+1,i) = norm(Lam_HP{l}(1,4:6) - Lam{i-1}(end ,4:6));
    end
  else
   dV2(1:MC_runs+1,i) = norm(Lam{i}(1,4:6) - Lam{i-1}(end,4:6));
   end
end

% final brake
dV3 = zeros(MC_runs+1, 1);
V_fin = zeros(1,3);  % Target velocity
if  disp_HP == np-1
    dV3(1,1) = norm( V_fin - Lam{np-1}(end,4:6));
    for l = 1: MC_runs
     Lam_HP{l} = sol_HP{disp_HP+l}.Lam_HP;
     dV3(l+1,1) = norm(V_fin - Lam_HP{l}(end,4:6));
    end
elseif disp_HP == 0
    dV3(1:MC_runs+1,1) = norm( V_fin - Lam{np-1}(end,4:6));
else
    dV3(1:MC_runs+1,1) = 0;
end
dV2_tot = sum(dV2,2);
dV_tot_adim = dV1 + dV2_tot + dV3 ;
dV_tot_dim = dV_tot_adim * cr3bp.L/cr3bp.T*2*pi; %km/s

% deltaV
output.dV.departure = dV1;
output.dV.HP = dV2;
output.dV.final = dV3;
output.dV.total_dim = dV_tot_dim;