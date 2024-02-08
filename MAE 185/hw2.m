clc
clear
load cylinder_re100
%% % w  = du/dx- dv/dy use first forward diff to find these quantities 
format long
timesteps = length(u);

%% find dx and dy
dx = mean(diff(x(:,1))); % average distance between x points
dy = mean(diff(y(1,:))); % y points

figure(1)
gif('vorticity_finitediff.gif')
for i = 1:timesteps
dvdx = ddx_central(squeeze(v(i,:,:)),dx);
dudy = ddy_central(squeeze(u(i,:,:)),dy);
w = dvdx-dudy;
pcolor(x,y,w);
rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]);
colorbar
title('$$\vec{w_z}$$', 'Interpreter', 'LaTeX')
xlabel('x')
ylabel('y')
axis equal tight
shading interp
drawnow
gif
end