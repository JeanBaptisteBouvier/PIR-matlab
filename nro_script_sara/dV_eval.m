function output = dV_eval(Lam, x0_rel, np, cr3bp, disp_HP, MC_runs, sol_HP)

% Initialization
% dV2 = zeros(np-1,1);
% dV3 = zeros(np-2,1);  %np-2 because before the first impulse from the orbit, the chaser has a non-zero velocity

% initial impulse from orbit
dV1 = norm(Lam{1}(1,4:6) - x0_rel(4:6)');
 
% check the presence of the dispersion
if disp_HP ==0
        k = 1:np-1;
        j = 2:np-1;
        MC_runs = 0;
else
        k = 1:disp_HP;
        j = 2:disp_HP;
end

% final brake for each hold point
dV2 = zeros(MC_runs+1, length(k));
for i= k
  if i == disp_HP
    dV2(1,i) =  norm(Lam{i}(end ,4:6));
    for l = 1: MC_runs
     Lam_HP{l} = sol_HP{i+l}.Lam_HP;
     dV2(l+1,i) = norm(Lam_HP{l}(end,4:6));
    end
  else
   dV2(1:MC_runs+1,i) =  norm(Lam{i}(end,4:6));

   end
end
 
% initial impulse from hold point
dV3 = zeros(MC_runs+1, length(k));
dV3(1:MC_runs+1,1) = dV1;
if np > 2
  for i= j
     if i == disp_HP 
         dV3(1,i) =  norm(Lam{i}(1 ,4:6));
         for l = 1: MC_runs
          Lam_HP{l} = sol_HP{i+l}.Lam_HP;
          dV3(l+1,i) =  norm(Lam_HP{l}(1 ,4:6));
         end
      else
          dV3(1:MC_runs+1,i) =  norm(Lam{i}(1,4:6));
      end
    end
else
    dV3 = 0;
end
 
dV_tot_adim = sum(dV2,2) + sum(dV3,2);
dV_tot_dim = dV_tot_adim * cr3bp.L/cr3bp.T*2*pi; %km/s

% deltaV
output.dV.departure = dV3;
output.dV.brake = dV2;
output.dV.total_dim = dV_tot_dim;