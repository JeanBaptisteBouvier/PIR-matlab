function [zrsyn] = rlvlh2rsyn_e(t, zlvlh, zrsyn_t, mu)
% Change of coordinates: from LVLH to RELATIVE Earth-Moon synodic.
%
% Inputs:
%  - zlvlh a 6*1 vector giving the LVLH state (with respect to the Earth).
%  - zrsyn_t a 6*1 vector giving the target state. 
%
% Outputs:
%  - zrlvlh a 6*1 vector giving the LVLH state (with respect to the Earth).
%

%--------------------------------------------------------------------------
% COC: from SYN to ECI for the Target
%--------------------------------------------------------------------------
zeci_t = syn2eci(t, zrsyn_t, mu);

%--------------------------------------------------------------------------
% COC: from LVLH to RCLI for the relative state
%--------------------------------------------------------------------------
zreci = rlvlh2reci(t, zlvlh, zeci_t, mu);

%--------------------------------------------------------------------------
% COC: from RLCI to RSYN for the relative state
%--------------------------------------------------------------------------
zrsyn = reci2rsyn(t, zreci);

end


