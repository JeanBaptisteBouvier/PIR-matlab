function [] = orbit_plot_lci( orbit, cr3bp, index )
% Plot an orbit in lunar-centered inertial coordinates

%--------------------------------------------------------------------------
%Create handle
%--------------------------------------------------------------------------
fig = figure(index);
hold on

%Default size
set(fig, 'defaultTextFontSize', 18);
set(fig, 'defaultAxesFontSize', 18);
set(fig, 'defaultTextFontWeight', 'bold');
set(fig, 'defaultTextHorizontalAlignment', 'center');
set(fig, 'defaultLineMarkerSize', 2);
set(fig, 'Position', [246 55 981 826]);

%--------------------------------------------------------------------------
%Settings
%--------------------------------------------------------------------------
axis equal
grid on
xlabel('X (-)')
ylabel('Y (-)')
zlabel('Z (-)')
title ('Trajectory in the lunar-centered inertial (LCI) frame');

%--------------------------------------------------------------------------
%Moon
%--------------------------------------------------------------------------
Rm2 = cr3bp.m2.Req/cr3bp.L;
[X_M3D, Y_M3D, Z_M3D] = sphere;
X_M3D = + Rm2*X_M3D;
Y_M3D = + Rm2*Y_M3D;
Z_M3D = + Rm2*Z_M3D;
% HMOON = surf(X_M3D, Y_M3D, Z_M3D, 'FaceColor', [23 153 179]./255, 'FaceLighting', 'none', 'EdgeColor', [9 63 87]./255);
HMOON = surf(X_M3D, Y_M3D, Z_M3D);

% Load Moon Image. CAREFUL: quite heavy for images!
load moonalb

% Set it on MOON
set(HMOON,'facecolor','texture',...
    'cdata',im2double(moonalb),...
    'edgecolor','none');
colormap(gray(256));


%--------------------------------------------------------------------------
%Plot
%--------------------------------------------------------------------------
plot3(orbit.ylci(:,1),orbit.ylci(:,2), orbit.ylci(:,3), 'Color',  rgb('dark green'), 'LineWidth', 1.5);

end

