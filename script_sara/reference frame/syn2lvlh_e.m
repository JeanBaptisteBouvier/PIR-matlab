function [zlvlh] = syn2lvlh_e(t, zsyn, zsyn_t, mu)
% Change of coordinates: from absolute Earth-Moon synodic LVLH.
%
% Inputs:
%  - zsyn a 6*1 vector giving the absolute state.
%  - zsyn_t a 6*1 vector giving the target state. 
%
% Outputs:
%  - zlvlh a 6*1 vector giving the LVLH state.
%

%--------------------------------------------------------------------------
% COC: from SYN to ECI for the Target
%--------------------------------------------------------------------------
zeci_t = syn2eci(t, zsyn_t, mu);

%--------------------------------------------------------------------------
% COC: from SYN to ECI for the absolute state
%--------------------------------------------------------------------------
zeci = syn2eci(t, zsyn, mu);

%--------------------------------------------------------------------------
% COC: from ECI to LVLH for the absolute state
%--------------------------------------------------------------------------
zlvlh = eci2lvlh(t, zeci, zeci_t, mu);

end


