%% 4.1
 
load('supersonicJetLES_xyPlane.mat');

cp = 1005 ; % J/kgK
cv = 718; %j/kgk
R = cp-cv;

U = prim2cons(rho,u,v,T,cv);

[rho2,u2,v2,T2,p2,e2,Et2] = cons2prim(U,R,cv);

mu_visoucous  = sutherland(T);

figure(1)

subplot(4,2,1)
pcolor(xx,yy,rho);
colorbar
shading interp
axis equal tight
title('\rho (kg/m^3)')
xlabel('x')
ylabel('y')

rho_u = squeeze(U(2,:,:)); 
subplot(4,2,2)
pcolor(xx,yy,rho_u);
colorbar
shading interp
axis equal tight
title('rho*u ({Kg}/{m^2s})',Interpreter='tex')
xlabel('x')
ylabel('y')

rho_v = squeeze(U(3,:,:)); 
subplot(4,2,3)
pcolor(xx,yy,rho_v);
colorbar
shading interp
axis equal tight
title('rho*v ({Kg}/{m^2s})',Interpreter='tex')
xlabel('x')
ylabel('y')

subplot(4,2,4)
pcolor(xx,yy,T);
colorbar
shading interp
axis equal tight
title('T (K)')
xlabel('x')
ylabel('y')

p = squeeze(p2);
subplot(4,2,5)
pcolor(xx,yy,p);
colorbar
shading interp
axis equal tight
title('P (Pa)')
xlabel('x')
ylabel('y')

e = squeeze(e2);
subplot(4,2,6)
pcolor(xx,yy,e);
colorbar
shading interp
axis equal tight
title('e ({J}/{kg})')
xlabel('x')
ylabel('y')

Et= squeeze(U(4,:,:)); 
subplot(4,2,7)
pcolor(xx,yy,Et);
colorbar
shading interp
axis equal tight
title('E_{t} ({J}/{m^3})',Interpreter='tex')
xlabel('x')
ylabel('y')

subplot(4,2,8)
pcolor(xx,yy,mu_visoucous);
colorbar
shading interp
axis equal tight
title('\mu ({Ns}/{m^2})',Interpreter='tex')
xlabel('x ')
ylabel('y')


function U = prim2cons(rho,u,v,T,cv)
%% Function prim2cons takes in 5 matrices (ngrid) representing primitive data about a test plane and returns a 4xmxn 
% where rho is density, u, v are x and y velocity

if(isequal(size(rho),size(u),size(v),size(T)))
    [nx,ny] =  size(rho) ;
    U = zeros(4,nx, ny); % create host vector U which has 4 elements, each a matrix of the same size as inputs
else
    sprintf("primitive matrices have incompatible sizes, all inputs must have equal dimensions");
end


rho_out = rho; % for some reason the function outputs rho too
U(1,:,:) = rho_out;

rho_u = rho.*u; 
U(2,:,:)= rho_u;

rho_v = rho.*v;
U(3,:,:) = rho_v;

e = cv*T;
Et = rho.*(e+(u.^2+v.^2)/2);
U(4,:,:)=Et;

end

function [rho,u,v,T,p,e,Et] = cons2prim(U,R,cv)
rho = squeeze(U(1,:,:));

rho_u  = squeeze(U(2,:,:));
u = rho_u./rho;

rho_v = squeeze(U(3,:,:));
v = rho_v./rho;

Et = squeeze(U(4,:,:));
e = Et./rho-(u.^2+v.^2)/2;
T = squeeze(e/cv); 

p = (rho*R).*T;
end
function [mu] = sutherland(T)
mu_0 = 1.735*10^(-5); % Ns/m^2

% assume t_0 == 300k? unindicated
T_0 = 300; 
S_1 = 110.4;

mu = mu_0*(T./T_0).^(3/2).*(T_0+S_1)./(T+S_1);


end

