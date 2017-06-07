function output = line_of_sight_corridor( alpha, beta, phi, theta, TOF , x0, min_dist)
% This Function compute the hold points in LVLH that are used in the Rendezvous code
% 
% alpha and beta are two trigger angles [rad]
% theta is the approach angle [rad]
% phi is an offset angle less than either trigger angle [rad]
% TOF is the time of flight [s]
% x0 is the position when the chaserfirst arrives at the corridor [km]
% min_dist is the minimum distance at which the chaser stop to maneuver [km]
A0 = norm(x0);
An= A0;
delta_X = abs(x0(1) - min_dist * cos((alpha+beta)/2)); % delta_X is the x-distance of the chaser from the target in LVLH that is the distance travelled in the corridor 
n=0;
hold_points = [];
delta_T = [];
while An > min_dist
    n = n+1;
%     An = A0 * (sin(phi)/sin(alpha+beta+phi))^n;
    An = An * (sin(phi)/sin(alpha+beta+phi));
    Bn = An * (sin(beta+alpha)/sin(phi));
    if mod(n,2)== 0
        hold_points = [hold_points; An*[cos(beta) 0 -sin(beta) ]];
        delta_T = [delta_T; Bn*abs(cos(alpha+phi))*TOF/delta_X];
    else 
        hold_points = [hold_points; An*[cos(alpha) 0 sin(alpha) ]];
        delta_T = [delta_T; abs(Bn*cos(beta+phi))*TOF/delta_X];
    end
end
output.hold_points =  hold_points;
output.delta_T = delta_T;
end