%% MIDTERM 2023
%Quinn, 
clear;
clc;
%% Generate nGrid, dx,dy
x_coords = linspace(0,1*10^(-5),75); % array containing each grid x coordinate
y_coords = linspace(0,8*10^(-6),80); % and y
[xx,yy]  = ndgrid(x_coords,y_coords); % generates x and y ndgrid coordinates for plotting
dx = mean(diff(xx(:,1))); % coul do analytically as well -\_(*_*)_/-
dy = mean(diff(yy(1,:)));
dt = 2.35*10^(-11); % timestep (s)
%% initialize initial Conditions - primitive variables
% datum labeled X_previous will store the value of primitive variable X at
% loop N in order to calulate X_next at loop number N+1

% Air Properties - constants
gamma = 1.4; 
R = 287 ; % J/Kg*K
T_0 = 288.1; %K
Pr = 0.71 ; % Prandtl No.
rho_0 = 1.225;%kg/m3
P_0 = 101300; % Pa
cv_0 = 718; % J/KgK
cp_0 = 1005: % J/kgK


% air properties - variable
T_previous = T_0*ones(75,80); % initial temperature grid (K)
a_previous = sqrt(gamma*R*T_previous);
rho_previous  = rho_0 *ones(75,80);
P_previous  = P_0*ones(75,80);
mu_previous = sutherland(T_last);
k_previous = 

u_previous = 4*a_previous; % X velocity is mach 4
v_previous = zeros(75,80); % no y velocity at start




%% Functions in Use

function [mu] = sutherland(T)
mu_0 = 1.735*10^(-5); % Ns/m^2
T_0 = 288.15; 
S_1 = 110.4;

mu = mu_0*(T./T_0).^(3/2).*(T_0+S_1)./(T+S_1);
end

function [rho,u,v,T,p,e,Et] = cons2prim(U,R,cv)
rho = U(1,:,:);

rho_u  = U(2,:,:);
u = rho_u./rho;

rho_v = U(3,:,:);
v = rho_v./rho;

Et = U(4,:,:);
e = Et./rho-(u.^2+v.^2)/2;
T = e/cv; 

p = (rho*R).*T;
end

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
