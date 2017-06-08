function output = line_of_sight_glide(alpha, phi, TOF, x0, min_dist)
% This Function compute the hold points in LVLH that are used in the Rendezvous code
% 
% alpha is the trigger angle [rad]
% phi is an offset angle less than either trigger angle [rad]
% TOF is the time of flight [s]
% x0 is the position when the chaserfirst arrives at the corridor [km]
% min_dist is the minimum distance at which the chaser stop to maneuver [km]
% 
% The calculation of the time needed between each hold point is a
% percentage of the TOF, corresponding to the ratio of the z-distance between
% the 2 hold points and the total z-distance
% 
% 

rel_angl = atan(x0(3)/x0(1));
if -alpha > rel_angl && rel_angl > alpha % outside
    error('The chaser is initially out of the corridor, redefine parameters') 
end



s = x0(1)/abs(x0(1)); % sign of x0(1)
hold_points = [s*abs(x0(3))/tan(alpha) 0 x0(3)]; % first hold point, chaser go straight, with the same z, until it reachs a wall of the corridor
A0 = norm(hold_points(1));

% to arrive from the good side
if x0(1)<0 && x0(3)<0
    A0 = -A0;
elseif x0(1)<0 && x0(3)>0
    A0 = -A0;
    alpha = -alpha;
    phi = - phi;
elseif x0(1)>0 && x0(3)<0
    alpha = -alpha;
    phi = - phi;
end

beta = alpha - phi;
An = A0;
D = abs(x0(3) - min_dist) + abs(x0(1)-hold_points(1,1)); % Distance travelled in the corridor 
n = 0;
delta_T = [abs(hold_points(1,1)-x0(1))*TOF/D];
while abs(An) > min_dist
    n = n+1;
    An = An * (sin(phi)/sin(alpha));
    Bn = An * (sin(beta)/sin(phi));
    hold_points = [hold_points; An*[cos(alpha + n*beta) 0 sin(alpha + n*beta) ]];
    delta_T = [delta_T; abs(hold_points(n+1,3)-hold_points(n,3))*TOF/D];
end
 

output.hold_points =  hold_points;
output.delta_T = delta_T;

figure(13)
hold on
grid on
plot([x0(1) hold_points(1,1)], [x0(3) hold_points(1,3)], 'b')
plot(x0(1), x0(3), '*g') % Chaser initial position
plot([A0*cos(alpha) -A0*cos(alpha)], [A0*sin(alpha) -A0*sin(alpha)], ':r') % Corridor
plot([A0*cos(alpha) -A0*cos(alpha)], [-A0*sin(alpha) A0*sin(alpha)], ':r') % Corridor
for k = 1:n-1
    plot(hold_points(k,1), hold_points(k,3), '*k')
    plot([hold_points(k,1) hold_points(k+1,1)], [hold_points(k,3) hold_points(k+1,3)], 'b')
end

end