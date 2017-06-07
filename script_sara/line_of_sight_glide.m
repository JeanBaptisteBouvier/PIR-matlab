function output = line_of_sight_glide(alpha, phi, theta, TOF, x0, min_dist)
% This Function compute the hold points in LVLH that are used in the Rendezvous code
% 
% alpha is the trigger angle [rad]
% theta is the approach angle [rad]
% phi is an offset angle less than either trigger angle [rad]
% TOF is the time of flight [s]
% x0 is the position when the chaserfirst arrives at the corridor [km]
% min_dist is the minimum distance at which the chaser stop to maneuver [km]

beta = alpha - phi; 
A0 = norm(x0);
An = A0;
delta_Z = abs(x0(3) - min_dist); % delta_X is the x-distance of the chaser from the target in LVLH that is the distance travelled in the corridor 
n=0;
hold_points = [];
delta_T = [];
while An > min_dist
    n = n+1;
    An = An * (sin(phi)/sin(alpha));
    Bn = An * (sin(beta)/sin(phi));
    hold_points = [hold_points; An*[cos(alpha + n*beta) 0 sin(alpha + n*beta) ]];
    delta_T = [delta_T; abs(Bn * sin(n*beta))*TOF/delta_Z];
 end
output.hold_points =  hold_points;
output.delta_T = delta_T;

end