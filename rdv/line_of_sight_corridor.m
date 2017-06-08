function output = line_of_sight_corridor( alpha, beta, phi, theta, TOF , x0, min_dist)
% This Function compute the hold points in LVLH that are used in the Rendezvous code
% 
% alpha and beta are two trigger angles [rad]
% theta is the approach angle [rad]
% phi is an offset angle less than either trigger angle [rad]
% TOF is the time of flight [s]
% x0 is the position when the chaserfirst arrives at the corridor [km]
% min_dist is the minimum distance at which the chaser stop to maneuver [km]


hold_points = [];
delta_T = [];
delta_X = abs(x0(1) - min_dist * cos((alpha+beta)/2)); % delta_X is the x-distance of the chaser from the target in LVLH that is the distance travelled in the corridor 

% Chech whether the Chaser is initially outside the corridor
rel_angl = atan(x0(3)/x0(1));
s = x0(1)/abs(x0(1)); % signe of x0(1)
if 0 < rel_angl && rel_angl < beta % in the upper corridor
    n=1;
    hold_points = [s*abs(x0(3))/tan(beta) 0 x0(3)];
    
elseif  0 > rel_angl && rel_angl > -alpha % in the lower corridor
    n=0;
    hold_points = [s*abs(x0(3))/tan(alpha) 0 x0(3)];
    
else % outside
    error('The chaser is initially out of the corridor, redefine parameters') 
end

A0 = norm(hold_points(1));
if hold_points(1,1)<0 % to arrive from the good side
    A0 = -A0;
end
An = A0;

while abs(An) > min_dist
    n = n+1;
    An = An * (sin(phi)/sin(alpha+beta+phi)); % recurrence
    Bn = An * (sin(beta+alpha)/sin(phi));
    if mod(n,2)== 0
        hold_points = [hold_points; An*[cos(beta) 0 -sin(beta) ]];
        delta_T = [delta_T; abs(Bn*cos(alpha+phi))*TOF/delta_X];
    else 
        hold_points = [hold_points; An*[cos(alpha) 0 sin(alpha) ]];
        delta_T = [delta_T; abs(Bn*cos(beta+phi))*TOF/delta_X];
    end
end

output.hold_points =  hold_points;
output.delta_T = delta_T;


figure
hold on
grid on
plot([A0*cos(alpha) -A0*cos(alpha)], [A0*sin(alpha) -A0*sin(alpha)], ':g')
plot([A0*cos(alpha) -A0*cos(alpha)], [-A0*sin(alpha) A0*sin(alpha)], ':g')
plot([x0(1) hold_points(1,1)], [x0(3) hold_points(1,3)], 'b')
for k = 1:n-1
    plot(hold_points(k,1), hold_points(k,3), '*k')
    plot([hold_points(k,1) hold_points(k+1,1)], [hold_points(k,3) hold_points(k+1,3)], 'b')
end


end