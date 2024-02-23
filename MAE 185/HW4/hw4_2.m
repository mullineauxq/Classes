%% HW 4.2

load('supersonicJetLES_xyPlane.mat');

cp = 1005 ; % J/kgK
cv = 718; %j/kgk
R = cp-cv;

dx = mean(diff(xx(:,1)));
dy  = mean(diff(yy(1,:)));

mu = sutherland(T); % create matrix of dynamic viscosities mu

dudx_fwd = ddx_fwd(u,dx); % compute first order 1st forward diff of u in x
dvdy_center = ddy_central(v,dy);% compute 2nd order 1st central diff in y 
divergence_velocity = dudx_fwd+dvdy_center; % divergence of velocity vector
dudy_central = ddy_central(u,dy); 
dvdx_fwd = ddx_fwd(v,dx);
Txx = 2*mu.*(dudx_fwd-1/3*divergence_velocity);% formula for normal shear Txx

Tyy = 2*mu.*(dvdy_center-1/3*divergence_velocity);% formula for normal shear Tyy

Txy = mu.*(dudy_central+dvdx_fwd); % fixed 


dudx_center = ddx_central(u,dx); % center diff of u in x 
dvdy_bwd = ddy_bwd(v,dy) ; % same procedure, with backwards diff of v in y
dudy_bwd = ddy_bwd(u,dy);
dvdx_central = ddx_central(v,dx);
divergence_velocity2 = dudx_center+dvdy_bwd; % divergence of velocity vector

Txx2 = 2*mu.*(dudx_center-1/3*divergence_velocity2);
Tyy2 = 2*mu.*(dvdy_bwd-1/3*divergence_velocity2 );

Txy2 = mu.*(dudy_bwd+dvdx_central);


%% Plotting 
figure(1)

subplot(2,3,1)
pcolor(xx,yy,Txx)
title('T_{xx}')
xlabel('x')
ylabel('y')
cb1 = colorbar;
ylabel(cb1,'Pa');
clim([-0.5,0.5]);
shading interp
axis equal tight 

subplot(2,3,2)
pcolor(xx,yy,Tyy)
title('T_{yy}')
xlabel('x')
ylabel('y')
cb1 = colorbar;
ylabel(cb1,'Pa');
clim([-0.5,0.5]);
shading interp
axis equal tight 

subplot(2,3,3)
pcolor(xx,yy,Txy)
title('T_{xy}')
xlabel('x');
ylabel('y');
cb1 = colorbar;
ylabel(cb1,'Pa');
clim([-0.5,0.5]);
shading interp
axis equal tight 

subplot(2,3,4)
pcolor(xx,yy,Txx2)
title('T_{xx}')
xlabel('x');
ylabel('y');
cb1 = colorbar;
ylabel(cb1,'Pa');
clim([-0.5,0.5]);
shading interp
axis equal tight 

subplot(2,3,5)
pcolor(xx,yy,Tyy2)
title('T_{yy}')
xlabel('x');
ylabel('y');
cb1 = colorbar;
ylabel(cb1,'Pa');
clim([-0.5,0.5]);
shading interp
axis equal tight 

subplot(2,3,6)
pcolor(xx,yy,Txy2)
title('T_{xy}')
xlabel('x');
ylabel('y');
cb1 = colorbar;
ylabel(cb1,'Pa');
clim([-0.5,0.5]);
shading interp
axis equal tight 


%% Functions
function [mu] = sutherland(T)
mu_0 = 1.735*10^(-5); % Ns/m^2

% assume t_0 == 300k? unindicated
T_0 = 300; 
S_1 = 110.4;

mu = mu_0*(T./T_0).^(3/2).*(T_0+S_1)./(T+S_1);


end