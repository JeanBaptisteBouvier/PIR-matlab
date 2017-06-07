% Michel Hénon : Numerical Exploration of the Restricted Problem V.
% List of 26 starting conditions of known DROs

Gamma = [6 5.5 5 4.5 4 3.5 3 2.5 2 1.5 1 0.5 0 -0.5 -1 -1.5 -2 -2.5 -3 -3.5 -4 -5 -6 -7 -8 -10];
Ksi = [-0.14779 -0.15888 -0.17169 -0.18661 -0.20421 -0.22523 -0.25071 -0.28212 -0.32163 -0.37252 -0.43991 -0.53182 -0.65966 -0.83185 -1.034 -1.23414 -1.41681 -1.58169 -1.73195 -1.87053 -1.99966 -2.23578 -2.44926 -2.64558 -2.82829 -3.16219];

% Transformation between two coordinates systems to obtain x0 and Co
% C0 not enough precise, I did more calculations to get the correction
mu = cr3bp.mu;
x0 = 1 - mu + ( mu^(1/3) ) * Ksi;
C0 = 3 + ( mu^(2/3) ) * Gamma - 4* mu - 2 * mu * Ksi.^3;
DRO = struct('x_0', x0,'C_j', C0);

load Dro_data;
load Dro_data_2; % is in reverse order, need to use fliplr()

figure;
hold on
grid on
C = horzcat(fliplr(Dro_data_2.C), Dro_data.C);
x = horzcat(fliplr(Dro_data_2.x), Dro_data.x);
T = horzcat(fliplr(Dro_data_2.T), Dro_data.T);

plot(x, C);
plot(x0, C0,'.')


% Draw the 26 known DRO contained in the DRO structur

% h = waitbar(0,'Computation in progress...');
% MaxIter = 26;
% for k = 1:MaxIter
%     if k==19
%         k = k+1;
%     end
%     waitbar(k/MaxIter);
%     x_0 = DRO.x_0(k);
%     Cj = DRO.C_j(k);
%     Ax = cr3bp.L * (1-cr3bp.mu - x_0);
%     dro = init_orbit(cr3bp, cr3bp.l2, cst.orbit.type.PLYAP, cst.orbit.family.PLANAR, Ax, cst);
%     v_y(k) = sqrt(abs( x_0^2  - Cj + 2*( (1-cr3bp.mu)/abs(x_0+cr3bp.mu) + cr3bp.mu/abs(x_0-1+cr3bp.mu) ) ));
%     yvg = [ x_0; 0; 0; 0 ; v_y(k) ; 0];
% 
%     orbit = orbit_refinement(cr3bp, dro, default, yvg, cst);
% end
% 
% close(h);
