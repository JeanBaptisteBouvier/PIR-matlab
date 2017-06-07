function out = bcfbp_derivatives_6(t, y, mu, theta0, ms, as, omegaS)
% BCFBP_DERIVATIVES_6 provide the equations of motion for the EARTH-MOON 
% BICIRCULAR FOUR-BODY PROBLEM (BCP) in MATLAB ODE format.
%
% OUT = BCFBP_DERIVATIVES_6(T, Y, MU, THETA0) computes the first-order
% equations of motion of the EARTH-MOON BCP at time T and state Y.
% The Earth-Moon mass ratio is MU. The initial position of the fourth body 
% (e.g. the Sun) is given by the scalar THETA0. The Sun mass, semi-major
% axis and mean motion must be provided by the user via the arguments MS,
% AS and OMEGAS.
%
% The equations of motion are available in equation (5.3.1) of Koon et al.
% "Dynamical Systems, the Three-Body Problem and Space Mission Design" 2006
% <a href="matlab: 
% web('http://www.cds.caltech.edu/~marsden/volume/missiondesign/KoLoMaRo_DMissionBook_2011-04-25.pdf','-browser')">(link)</a>. 
%
% See also BCFBP_DERIVATIVES_42
%
% BLB 2015

%Output declaration
out = (1:6)';

%--------------------------------------------------------------------------
% Update first derivatives of the potential \bar{U} (cf Koon et al. 2006)
%--------------------------------------------------------------------------
d1_ub = d1_u_b4(mu, ms, as, y(1),y(2),y(3), theta0+omegaS*t);

%--------------------------------------------------------------------------
%Phase space derivatives
%--------------------------------------------------------------------------
out(1) = y(4);
out(2) = y(5);
out(3) = y(6);
out(4) = d1_ub(1) + 2*y(5);
out(5) = d1_ub(2) - 2*y(4);
out(6) = d1_ub(3);

end


%--------------------------------------------------------------------------
% Subroutines
%--------------------------------------------------------------------------
function xout = d1_u_b4(mu, ms, as, x,y,z, theta)
% First derivative of the energy potential
%
% BLB 2015

r1 = sqrt((x+mu)^2 + y^2 + z^2);
r2 = sqrt((x-1+mu)^2 + y^2 + z^2);
rs = sqrt((x-as*cos(theta))^2 + (y-as*sin(theta))^2 + z^2);
xout(1) = x - (1-mu)/r1^3*(x+mu) - mu/r2^3*(x-1+mu) - ms/rs^3*(x-as*cos(theta)) - ms/as^2*cos(theta);
xout(2) = y - (1-mu)/r1^3*y      - mu/r2^3*y        - ms/rs^3*(y-as*sin(theta)) - ms/as^2*sin(theta);
xout(3) =   - (1-mu)/r1^3*z      - mu/r2^3*z        - ms/rs^3*z;
end