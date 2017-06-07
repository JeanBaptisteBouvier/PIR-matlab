% This algorithme choose the best targeting algorithm between CW, SL and LR, 
% according to our data base of errors. It has been computed with a gap of
% 2000 km between the chaser and the target, and an angular gap of 10°
% 
% Two parameters can be chosen : the Ax_Target and the Theta_Target
% 
% When the Ax and the angle requested is not in our data base,
% interpolations are lead to find the result. This algorithm  compute a 2D
% or a 3D interpolation wether one or two of the parameters does not belong
% to our data.
% 
% 
% 
% 
% 
% 
% By JBB

init;
cr3bp = init_CR3BP('EARTH', 'MOON', default);
load Comp

%% Parameters

Ax_Target = 180100; %radius of the orbit, between 16'000km and 250'000km
Theta_Target = 165; %in degrees


%% Computation

size_inf = 16000; % radius of the smallest target orbit (km)
size_step = 1000; % step between each target orbit (km)
size_sup = 250000; % radius of the biggest target orbit (km)

angl_inf = 0; % smallest angular position of the chaser (degree)
angl_step = 2; % step between each angles (degree)
angl_sup = 345; % biggest angular position of the chaser (degree)

Ax = size_inf:size_step:size_sup-size_step ;
Theta = angl_inf:angl_step:angl_sup ;
[T,S] = meshgrid(Theta, Ax);

% Interpolation
intp_size_stp = gcd(mod(Ax_Target, size_step), size_step); % the step of interpolation can be smaller than the one of calculus
intp_angl_stp = gcd(mod(Theta_Target, angl_step),angl_step);

% Construction of the list Err = [CW_error, SL_error, LR_error]

if intp_size_stp > 0 && intp_angl_stp > 0 % both steps are too loose, a 3D interpolation is needed
    
    Angl = angl_inf:intp_angl_stp:angl_sup;
    Size = size_inf:intp_size_stp:size_sup;
    [Angl, Size] = meshgrid(Angl, Size);
    L_CW = interp2(S, T, Comp.CW, Size, Angl, 'spline');
    L_SL = interp2(S, T, Comp.SL, Size, Angl, 'spline');
    L_LR = interp2(S, T, Comp.LR, Size, Angl, 'spline');
    Err(1) = L_CW(Ax_Target/intp_size_stp, Theta_Target/intp_angl_stp); %interpolation and retrieve the one value we are intersted in
    Err(2) = L_SL(Ax_Target/intp_size_stp, Theta_Target/intp_angl_stp);
    Err(3) = L_LR(Ax_Target/intp_size_stp, Theta_Target/intp_angl_stp);
    
    
elseif intp_size_stp == 0 && intp_angl_stp > 0 % only the angle step is too loose
    
    size_index = (Ax_Target-size_inf)/size_step +1;
    Angl = angl_inf:intp_angl_stp:angl_sup;
    L_CW = interp1(Theta, Comp.CW(size_index,:), Angl, 'spline');
    L_SL = interp1(Theta, Comp.SL(size_index,:), Angl, 'spline');
    L_LR = interp1(Theta, Comp.LR(size_index,:), Angl, 'spline');
    Err(1) = L_CW(Theta_Target/intp_angl_stp);
    Err(2) = L_SL(Theta_Target/intp_angl_stp);
    Err(3) = L_LR(Theta_Target/intp_angl_stp);
    
elseif intp_size_stp > 0 && intp_angl_stp == 0 % only the size step is too loose
    
    angl_index = (Theta_Target-angl_inf)/angl_step +1;
    Size = size_inf:intp_size_stp:size_sup;
    L_CW = interp1(Ax, Comp.CW(:,angl_index), Size, 'spline');
    L_SL = interp1(Ax, Comp.SL(:,angl_index), Size, 'spline');
    L_LR = interp1(Ax, Comp.LR(:,angl_index), Size, 'spline');
    Err(1) = L_CW(Ax_Target/intp_size_stp);
    Err(2) = L_SL(Ax_Target/intp_size_stp);
    Err(3) = L_LR(Ax_Target/intp_size_stp);
    
else % no interpolation is needed, the value is already in our database
    
    size_index = (Ax_Target-size_inf)/size_step +1;
    angl_index = (Theta_Target-angl_inf)/angl_step +1;
    Err(1) = Comp.CW(size_index,angl_index);
    Err(2) = Comp.SL(size_index,angl_index);
    Err(3) = Comp.LR(size_index, angl_index);

end



Err = Err*cr3bp.L;  % Conversion in kilometers

if Err(3)<Err(2) && Err(3)<Err(1)
    sprintf(['The LR algo is the best, with an error of ' num2str(Err(3)) ' km'])
elseif Err(2)<Err(3)
    sprintf(['The SL algo is the best, with an error of ' num2str(Err(2)) ' km'])
else
    sprintf(['The CW algo is the best, with an error of ' num2str(Err(1)) ' km'])
end

