function [zrlvlh] = rsyn2rlvlh_e(t, zrsyn, zrsyn_t, mu)
% Change of coordinates: from RELATIVE Earth-Moon synodic LVLH (with respect to the Earth).
%
% Inputs:
%  - zrsyn a 6*1 vector giving the relative state.
%  - zrsyn_t a 6*1 vector giving the target state. 
%
% Outputs:
%  - zrlvlh a 6*1 vector giving the LVLH state.
%

%--------------------------------------------------------------------------
% COC: from SYN to ECI for the Target
%--------------------------------------------------------------------------
zeci_t = syn2eci(t, zrsyn_t, mu);

%--------------------------------------------------------------------------
% COC: from RSYN to RECI for the relative state
%--------------------------------------------------------------------------
zreci = rsyn2reci(t, zrsyn);

%--------------------------------------------------------------------------
% COC: from RECI to LVLH for the relative state
%--------------------------------------------------------------------------
zrlvlh = reci2rlvlh(t, zreci, zeci_t, mu);

end


