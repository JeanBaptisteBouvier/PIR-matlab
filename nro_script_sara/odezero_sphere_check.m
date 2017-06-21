function [value,isterminal,direction] = odezero_sphere_check(t,y,r_sphere)
%
% [VALUE,ISTERMINAL,DIRECTION] = ODEZERO_SPHERE_CHECK(T, Y, EVENT) 
% is an event routine in MATLAB ODE format (see event in MATLAB help). 
% The condition of the event triggering is:
%     sqrt(x_pos(1)^2+x_pos(2)^2+x_pos(3)^2) = r_sphere;
%
% y [km] is the chaser position at the t [s] instant 

%Event parameters

value = sqrt(y(1)^2+ y(2)^2+ y(3)^2) - r_sphere;
isterminal = 0;
direction = -1;

end