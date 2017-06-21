function [zlvlh] = reci2rlvlh(t, zreci, zeci_t, mu)
% Change of coordinates: from relative ECI to LVLH  coordinates.

%--------------------------------------------------------------------------
% Rotation matrix
%--------------------------------------------------------------------------
Cs = lvlh2eciRotMat(zeci_t, t, mu);

%--------------------------------------------------------------------------
% COC
%--------------------------------------------------------------------------
zlvlh = Cs\zreci;
end
