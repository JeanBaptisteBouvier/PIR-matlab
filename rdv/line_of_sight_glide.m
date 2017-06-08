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

rel_angl = atan(x0(3)/x0(1));
if -alpha > rel_angl && rel_angl > alpha % outside
    error('The chaser is initially out of the corridor, redefine parameters') 
end

s = x0(1)/abs(x0(1)); % signe of x0(1)
hold_points = [s*abs(x0(3))/tan(alpha) 0 x0(3)]; % 
A0 = norm(hold_points(1));

% to arrive from the good side
if x0(1)<0 && x0(3)<0
    A0 = -A0;
    
elseif x0(1)<0 && x0(3)>0
    A0 = -A0;
    alpha = -alpha;
elseif x0(1)>0 && x0(3)<0
    alpha = -alpha;
else
    
end
An = A0;
delta_Z = abs(x0(3) - min_dist); % delta_X is the x-distance of the chaser from the target in LVLH that is the distance travelled in the corridor 
n=0;
hold_points = [];
delta_T = [];
hold_points = [hold_points; An*[cos(alpha + n*beta) 0 sin(alpha + n*beta) ]];
% hold_points = [hold_points; An*[1 0 1]./sqrt(2) ];
while abs(An) > min_dist
    n = n+1;
    An = An * (sin(phi)/sin(alpha));
    Bn = An * (sin(beta)/sin(phi));
    hold_points = [hold_points; An*[cos(alpha + n*beta) 0 sin(alpha + n*beta) ]];
    delta_T = [delta_T; abs(Bn * sin(n*beta))*TOF/delta_Z];
end
 

output.hold_points =  hold_points;
output.delta_T = delta_T;

figure
hold on
grid on
plot([x0(1) hold_points(1,1)], [x0(3) hold_points(1,3)], 'b')
plot(x0(1), x0(3), '*g')
for k = 1:n-1
    plot(hold_points(k,1), hold_points(k,3), '*k')
    plot([hold_points(k,1) hold_points(k+1,1)], [hold_points(k,3) hold_points(k+1,3)], 'b')
end

end