function output = sphere_plot(r, x0, y0,z0)
%x0, y0 and z0 are the coordinates of the center of the sphere
%r is the radius of the sphere
[x,y,z] = sphere(50);
x = x*r + x0;
y = y*r + y0;
z = z*r + z0;
lightGrey = 0.8*[1 1 1]; % It looks better if the lines are lighter
output = surface(x,y,z,'FaceColor', 'none','EdgeColor',lightGrey)
end
