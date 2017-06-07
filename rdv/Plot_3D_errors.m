% Plot the 3D comparison between the 3 algorithms on the radius of the
% orbits and the angular position of the Chaser
% 
% launch Targeting_Comparison.m to create the Comp structur
init;
cr3bp = init_CR3BP('EARTH', 'MOON', default);
load Comp

gap = 2000; % gap between the chaser and target orbits (km)
size_inf = 16000; % radius of the smallest target orbit (km)
size_step = 1000; % step between each target orbit (km)
size_sup = 250000; % radius of the biggest target orbit (km)

angl_inf = 0; % smallest angular position of the chaser (degree)
angl_step = 2.5; % step between each angles (degree)
angl_sup = 345; % biggest angular position of the chaser (degree)

Ax = size_inf:size_step:size_sup-size_step ;
Theta = angl_inf:angl_step:angl_sup ;
Theta_cut = horzcat(angl_inf:angl_step:165, 182.5:angl_step:angl_sup);

% Use of interpolation to find the points around 180°
l = length(Comp.CW(:,1));
L = length(Theta);
for k = 1:l
    CW(k,:) = interp1(Theta_cut, Comp.CW(k,:)*cr3bp.L, Theta, 'spline');
    SL(k,:) = interp1(Theta_cut, Comp.SL(k,:)*cr3bp.L, Theta, 'spline');
    LR(k,:) = interp1(Theta_cut, Comp.LR(k,:)*cr3bp.L, Theta, 'spline');
end

figure(1);
hold on
grid on
[Theta,Size] = meshgrid(Theta, Ax);

surf(Size,Theta,CW, ones(l,L));
surf(Size,Theta,SL, zeros(l,L));
surf(Size,Theta,LR, 0.5*ones(l,L)); 

legend('Clohessy-Wiltshire', 'Straight-Line', 'Linearized Relative')
title('Error comparison between Targeting algorithms')
xlabel('Radius of the Target DRO, which is  2 000 km less than the radius of the Chaser DRO')
ylabel('Initial angular position of the Chaser, 10 degrees behind the Target')
zlabel('Error in km')
hold off