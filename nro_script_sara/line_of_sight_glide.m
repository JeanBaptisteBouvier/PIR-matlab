function output = line_of_sight_glide(alpha, beta, phi, TOF, x0, min_dist, cr3bp)
% This Function compute the hold points in LVLH that are used in the Rendezvous code
% 
% alpha is the trigger angle [rad]
% phi is an offset angle less than either trigger angle [rad]
% TOF is the time of flight [s]
% x0 is the position when the chaserfirst arrives at the corridor [km]
% min_dist is the minimum distance at which the chaser stop to maneuver [km]

 
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
delta_Z = abs(x0(3) - min_dist); % delta_Z is the z-distance of the chaser from the target in LVLH that is the distance travelled in the corridor 
n = 0;
delta_T = [abs(hold_points(1,3))*TOF/delta_Z];
while abs(An) > min_dist
    n = n+1;
    An = An * (sin(phi)/sin(alpha));
    Bn = An * (sin(beta)/sin(phi));
    hold_points = [hold_points; An*[cos(alpha + n*beta) 0 sin(alpha + n*beta) ]];
    delta_T = [delta_T; abs(Bn * sin(n*beta))*TOF/delta_Z];
end
 

output.hold_points =  hold_points;
output.delta_T = delta_T;

figure(13)
hold on
grid on
p1 = plot([cr3bp.L*x0(1) cr3bp.L*hold_points(1,1)], [cr3bp.L*x0(3) cr3bp.L*hold_points(1,3)], 'b','LineWidth',1);
p2 = plot(cr3bp.L*x0(1), cr3bp.L*x0(3), 'go', 'LineWidth',1); % Chaser initial position
p3 = plot([cr3bp.L*A0*cos(alpha) -cr3bp.L*A0*cos(alpha)], [cr3bp.L*A0*sin(alpha) -cr3bp.L*A0*sin(alpha)], ':c','LineWidth',1.5); % Segment associated to the trigger angle
p4 = plot([cr3bp.L*A0*cos(alpha) -cr3bp.L*A0*cos(alpha)], [-cr3bp.L*A0*sin(alpha) cr3bp.L*A0*sin(alpha)], ':c','LineWidth',1.5); % Segment associated to the trigger angle
p7 = plot(0,0,'ko','LineWidth',1);
for k = 1:n+1
    p5(k) = plot(cr3bp.L*hold_points(k,1), cr3bp.L*hold_points(k,3), 'mo','LineWidth',1);
    str{k} = cellstr(sprintf('%s %d','HP', k));
end
for k = 1:n
    p6 = plot([cr3bp.L*hold_points(k,1) cr3bp.L*hold_points(k+1,1)], [cr3bp.L*hold_points(k,3) cr3bp.L*hold_points(k+1,3)], 'b', 'LineWidth',1);
end
  pp = [p2 p3 p5 p6 p7];
  nom = ['Chaser initial position', 'Segment associated to the trigger angle', str{:} ,'Trajectory' ,'Docking'];
  legend(pp,nom)
  title('Rendezvous Trajectory Design in the LVLH Frame')
  xlabel('Downrange [km]')
  ylabel('Altitude [km]')
  
end