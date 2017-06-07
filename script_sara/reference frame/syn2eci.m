function [zeci] = syn2eci(t, zsyn, mu)
% Change of coordinates: from Earth-Moon synodic to Earth-centered inertial
% coordinates.
%--------------------------------------------------------------------------
% Earth state in EM synodical coordinates
%--------------------------------------------------------------------------
zsyn_e = [-mu ; 0 ; 0 ;  0 ; 0 ; 0 ];

%--------------------------------------------------------------------------
% COC
%--------------------------------------------------------------------------
zeci = eciRotMat(t)*(zsyn + zsyn_e);

end


%--------------------------------------------------------------------------
% Subroutines
%--------------------------------------------------------------------------
function RotMat = eciRotMat(theta)
% Compute the rotation matrix associated to the angle theta
%
%   R =  | R11    0  |
%        | R21  R11  |
% with
%
%         | c -s 0 |          | -s -c 0 |
%   R11 = | s  c 0 |,   R21 = |  c -s 0 | 
%         | 0  0 1 |          |  0  0 0 |
% and
%       c = cos(theta), s = sin(theta)
%
% BLB 2016
RotMat11 = [+cos(theta) -sin(theta) 0; +sin(theta) +cos(theta) 0 ; 0 0 1];
RotMat21 = [-sin(theta) -cos(theta) 0; +cos(theta) -sin(theta) 0 ; 0 0 0];
RotMat   = [RotMat11 zeros(3) ; RotMat21 RotMat11];
end

