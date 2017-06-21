function [ time, transfer, t_events_AS, y_events_AS, t_events_KOS, y_events_KOS] = Dispersion(V_plus, V_minus, x_HP, t_HP, x_T,  t_disp , err_mod, err_angle, phi0, cr3bp, options, options2)
% This function wil create the dispersion vector of the velocity relative
% to an hold point and it will propagate this dispersion along the
% rendezvous trajectory
%
% V_plus is the velocity that the chaser has when it arrives in the HP [km/s]
% V_plus is the velocity that the chaser has when it starts from the HP [km/s]
% x_HP is the position vector of the HP [km]
% t_disp represents the time of propagation [s]
% t_HP represents the time that the chaser needs to arrive at the selected HP [s]
% x_T is the vector relative to position and velocity of the target
% err_mod represents the percentage of the execution error (3sigma) in magnitude
% err_angle represents the execution error (3sigma) in pointing [°]
% phi0 is the initial state transition matrix for the Lambert problem

delta_V = V_plus- V_minus; % velocity vector relative to the maneuver in the HP in which the dispersion will be added
mean_mod = norm(delta_V);  % mean of the velocity magnitude
n = delta_V/norm(delta_V); % direction of the velocity vector

% Dispersion generation of the velocity magnitude
sigma_mod = err_mod/3*mean_mod;
SNR_mod = mean_mod/sigma_mod;
delta_V_mod_new = awgn(mean_mod, SNR_mod); 

% Dispersion generation of the velocity direction 
eps = err_angle^2/(2*(n(1)+ n(2)+ n(3)));
sigma_angle= abs(eps)/3;
theta = err_angle + 1; % initial random guess to enter in the while-cycle
while theta > err_angle  %[°]
    disper_angle_1 = normrnd(0,sigma_angle); 
    disper_angle_2 = normrnd(0,sigma_angle);  
    disper_angle_3 = normrnd(0,sigma_angle);  
    n_new = [n(1) + disper_angle_1, n(2) + disper_angle_2, n(3) + disper_angle_3];
    theta = atan2(norm(cross(n,n_new)), dot(n,n_new))*180/pi;  % angle [°] between the old and the new direction vector
end
% Dispersion propagation
delta_V_new = delta_V_mod_new.*n_new;
V_new = delta_V_new + V_minus;
x0 = [x_HP, V_new];
t_disp_norm = t_disp*2*pi/cr3bp.T;
 
[time, transfer, t_events_AS, y_events_AS, ~] = ode113(@(t,y)crtbp(t,y,cr3bp.mu),[t_HP  t_disp_norm+t_HP],[x0, x_T, phi0']', options); 

% check dispersion KOS
r_sphere = 2/cr3bp.L;   % radius of the Approach Sphere to ensure trajectory safety policies [km]
for i = 1:size(transfer,1);
    radius(i) = sqrt(transfer(i,1)^2+ transfer(i,2)^2+ transfer(i,3)^2)
    if radius(i) <= r_sphere % verify if we are inside oroutside the AS in order to check the condition on KOS
       [time, transfer, t_events_KOS, y_events_KOS, ~] = ode113(@(t,y)crtbp(t,y,cr3bp.mu),[t_HP  t_disp_norm+t_HP],[x0, x_T, phi0']', options2);
    else
       t_events_KOS = [];
       y_events_KOS = [];
    end
end

end
