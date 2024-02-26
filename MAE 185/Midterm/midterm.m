%% MIDTERM 2023
%Quinn, 

%% WHAT NEEDS FIXING? ADDS NOTES HERE
% prompt says to avoid calculating the value of all primitives at
% boundaries. This code calculates at boundaries unfortunately. 

% figure out why all primitives have a zero imaginary component. its
% messing with plotting
%%
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
cv = 718; % J/KgK
cp = 1005; % J/kgK


% air properties - variable
T_previous = T_0*ones(75,80); % initial temperature grid (K)
a_previous = real(sqrt(gamma*R*T_previous)); % local speed of sound at each point
rho_previous  = rho_0 *ones(75,80); % density
P_previous  = P_0*ones(75,80); % pressure
mu_previous = sutherland(T_previous); % sutherland's law viscosity
k_previous = cp/Pr*mu_previous; % conductivity
e_previous = cv*T_previous; % specific internal energy

u_0 = 4*mean(mean(a_previous)); % set initial x velocity to mach 4 (used for setting BCs later)

u_previous = 4*a_previous; % X velocity is mach 4
v_previous = zeros(75,80); % no y velocity at start
Et_previous = rho_previous.*(e_previous+1/2*(u_previous.^2+v_previous.^2));

% Create conservative vectors/ preallocate others 
U_previous = prim2cons(rho_previous, u_previous,v_previous,T_previous, cv); % Stores U, the conservative elements vector
NS_size = size(U_previous); % 
E_previous = zeros(NS_size);
F_previous = zeros(NS_size);
U_bar = zeros(NS_size);
E_bar = zeros(NS_size);
F_bar = zeros(NS_size);
U_next = zeros(NS_size); % I may go back and make this just U and remove the whole '_previous' shit. Doesnt seem necessary
dEdx_fwd = zeros(NS_size);% preallocating seems to avoid some problems
dFdy_fwd = zeros(NS_size);
dEbardx_bwd = zeros(NS_size);
dFbardy_bwd = zeros(NS_size);

%% Time Loop 

%%gif('developing_flow_M4.gif')
for n = 1:1500
   %% begin each loop by updating the vectors E and F with the values provided by U_previous
    % To do so, collect primitive variables,assuming all BCs have been
    % enforced.
    [rho_previous,u_previous,v_previous,T_previous,P_previous,e_previous,Et_previous] = cons2prim(U_previous,R,cv);
%%
    % 1: Calculate Derivatives FOR E (Predictor) (BACKWARDS IN X CENTRAL IN Y)

    dudx_bwd  = ddx_bwd(u_previous,dx); % first backward diff of u in x
    dvdx_bwd = ddx_bwd(v_previous,dx);
    dTdx_bwd = ddx_bwd(T_previous,dx);
    dvdy_central = ddy_central(v_previous,dy);
    dudy_central  = ddy_central(u_previous,dy);

    %1.1 - calculate shears and heat flux
    Tau_xx = 2*mu_previous.*(dudx_bwd-1/3*(dudx_bwd+dvdy_central));
    Tau_xy = mu_previous.*(dudy_central+dvdx_bwd);
    qx  = -k_previous.*dTdx_bwd; 

    % 2: Create E, calculate Derivative OF E (Forward in X)
    
    E_previous(1,:,:) = U_previous(2,:,:);
    E_previous(2,:,:) = rho_previous.*u_previous.^2 +P_previous-Tau_xx;
    E_previous(3,:,:) = rho_previous.*u_previous.*v_previous -Tau_xy;
    E_previous(4,:,:) = (Et_previous+P_previous).*u_previous-u_previous.*Tau_xx-v_previous.*Tau_xy+qx; 

   
    % to calculate dEdx_fwd I used a loop through each of the 4
    % conservative E fields. Seems to avoid some coding problems
    
    for j = 1:4
       dEdx_fwd(j,:,:) = ddx_fwd(squeeze(E_previous(j,:,:)),dx);
    end

    % 3: Calculate Derivatives FOR F (predictor)(BACKWARDS IN Y CENTRAL IN X)
    dudy_bwd = ddy_bwd(u_previous,dy);
    dvdy_bwd = ddy_bwd(v_previous,dy);
    dTdy_bwd = ddy_bwd(T_previous,dy); 

    dvdx_central  = ddx_central(v_previous,dx);
    dudx_central = ddx_central(u_previous,dx);
 
    % 3.1 Calculate shears and heat flux in y
    Tau_xy  =mu_previous.*(dudy_bwd+dvdx_central); % overwrite Tau_xy to change bias
    Tau_yy = 2*mu_previous.*(dvdy_bwd -1/3*(dudx_central+dvdy_bwd));
    qy = -k_previous.*dTdy_bwd;

    % 4: Calculate Derivative OF F (Forward in Y)

    F_previous(1,:,:) = U_previous(3,:,:);
    F_previous(2,:,:) = rho_previous.*u_previous.*v_previous - Tau_xy;
    F_previous(3,:,:) = rho_previous.*v_previous.^2+P_previous-Tau_yy;
    F_previous(4,:,:) = (Et_previous+P_previous).*v_previous-u_previous.*Tau_xy-v_previous.*Tau_yy+qy;
   
    % Just as with E, loop thru each of F's 4 fields to find the forward
    % deriv in x. avoids some errors
    for j = 1:4
       dFdy_fwd(j,:,:) = ddy_fwd(squeeze(F_previous(j,:,:)),dy);
    end

    % 5: Calculate Ubar Using derivative of E and F using vector column form
    % of Navier Stokes

    U_bar = U_previous-dt*(dEdx_fwd)-dt*(dFdy_fwd); % fwd diffs in predictor step
    
    %5.1 re calculate all primitive variables then enforce BCs
    [rho_previous,u_previous,v_previous,T_previous,P_previous,e_previous,Et_previous] = cons2prim(U_bar,R,cv); %% just overwriting the timestep n vars. wont be needed again
    
    % 5.1.1 BCs - WALL & leading edge  (First may be overwritten by leading edge BC)
    T_previous(:,1) = T_0;
    P_previous(:,1) = extrapolate2(P_previous,'bottom'); % 2nd order extrap. @ wall
    P_previous(1,1)= P_0;
    u_previous(:,1) = 0; % apply no-slip
    u_previous(1,1) = u_0;
    v_previous(:,1) = 0; % solid wall condition
    rho_previous(:,1) = P_previous(:,1)./(R*T_previous(:,1)); % set rho BCs using ideal gas law
    rho_previous(1,1) = rho_0 ; % This BC is not explicitly stated, but we must assume a BC for rho if P and T are fxed and P = rho * R *T
    e_previous(:,1) = cv*T_previous(:,1); % set sp. internal energy bcs using calorically perfect assumption
    Et_previous(:,1) = rho_previous(:,1).*(e_previous(:,1)+1/2*u_previous(:,1).^2);
    mu_previous(:,1) = sutherland(T_previous(:,1)); % recalculate viscosity with Temp BC in mind. 
    k_previous(:,1) = cp/Pr*mu_previous(:,1); % recalculate therm. conductivity with BC values of mu in mind

    %5.1.2 BCs - INLET Boundary (Will not Overwrite Wall)
    T_previous(1,2:end) = T_0;
    P_previous(1,2:end) = P_0;
    u_previous(1,2:end) = u_0;
    v_previous(1,2:end) = 0 ;
    rho_previous(1,2:end) = P_0/(R*T_0);
    e_previous(1,2:end) = cv*T_0;
    Et_previous(1,2:end) = rho_previous(1,2:end).*(e_previous(1,2:end)+1/2*(u_previous(1,2:end).^2+v_previous(1,2:end).^2));
    mu_previous(1,2:end) = sutherland(T_previous(1,2:end));
    k_previous(1,2:end) = cp/Pr*mu_previous(1,2:end);

    %5.1.3 - Far Field BCs (same condition as inlet, but on top boundary)
    %(overwrites outlet BC at corner)

    T_previous(2:end,end) = T_0;
    P_previous(2:end,end) = P_0;
    u_previous(2:end,end) = u_0;
    v_previous(2:end,end) = 0 ;
    rho_previous(2:end,end) = P_0/(R*T_0);
    e_previous(2:end,end) = cv*T_0;
    Et_previous(2:end,end) = rho_previous(2:end,end).*(e_previous(2:end,end)+1/2*(u_previous(2:end,end).^2+v_previous(2:end,end).^2));
    mu_previous(2:end,end) = sutherland(T_previous(2:end,end));
    k_previous(2:end,end) = cp/Pr*mu_previous(2:end,end);

    % 5.1.4 Outlet BCS (as specified by prompt) will not overwrite far
    % field or wall conditions!
    
    T_previous(end,2:(end-1)) = extrapolate2(T_previous,'right');
    P_previous(end,2:(end-1)) = extrapolate2(P_previous,'right');
    u_previous(end,2:(end-1)) = extrapolate2(u_previous,'right');
    v_previous(end,2:(end-1)) = extrapolate2(v_previous,'right');
    rho_previous(end,2:(end-1)) = P_previous(end,2:(end-1))./(R*T_previous(end,2:(end-1)));
    e_previous(end,2:(end-1)) = cv*T_previous(end,2:(end-1));
    Et_previous(end,2:(end-1)) = rho_previous(end,2:(end-1)).*(e_previous(end,2:(end-1))+1/2*(u_previous(end,2:(end-1)).^2+v_previous(end,2:(end-1)).^2));
    mu_previous(end,2:(end-1)) = sutherland(T_previous(end,2:(end-1)));
    k_previous(end,2:(end-1)) = cp/Pr*mu_previous(end,2:(end-1));

    % 5.2 update Ubar with boundary conditions
    U_bar = prim2cons(rho_previous,u_previous,v_previous,T_previous,cv);
    %At this stage, alll BCs are enforced in the predictor step, and Ubar
    %is updated. E bar and F bar must be recalculated

    % 6: Calculate Derivatives FOR Ebar (forwards in X, Central in y)
    dudx_fwd  = ddx_fwd(u_previous,dx); % first backward diff of u (with BCs updated) in x
    dvdx_fwd = ddx_fwd(v_previous,dx);
    dTdx_fwd = ddx_fwd(T_previous,dx);
    dvdy_central = ddy_central(v_previous,dy);
    dudy_central  = ddy_central(u_previous,dy);

    %6.1 - calculate shears and heat flux for E
    Tau_xx = 2*mu_previous.*(dudx_fwd-1/3*(dudx_fwd+dvdy_central));
    Tau_xy = mu_previous.*(dudy_central+dvdx_fwd);
    qx  = -k_previous.*dTdx_fwd; 
    
    % 7 - Create E_bar and take its derivative (backwards in x)
    E_bar(1,:,:) = U_bar(2,:,:);
    E_bar(2,:,:) = rho_previous.*u_previous.^2 +P_previous-Tau_xx;
    E_bar(3,:,:) = rho_previous.*u_previous.*v_previous -Tau_xy;
    E_bar(4,:,:) = (Et_previous+P_previous).*u_previous-u_previous.*Tau_xx-v_previous.*Tau_xy+qx; 
    
    for j = 1:4
       dEbardx_bwd(j,:,:) = ddx_bwd(squeeze(E_bar(j,:,:)),dx);
    end

    % 8: Calculate  Derivatives FOR Fbar (forwards in y, central in x)
    dudy_fwd = ddy_fwd(u_previous,dy);
    dvdy_fwd = ddy_fwd(v_previous,dy);
    dTdy_fwd = ddy_fwd(T_previous,dy); 
    
    dvdx_central  = ddx_central(v_previous,dx);
    dudx_central = ddx_central(u_previous,dx);

    % 8.1 Calculate Shears for Fbar
    Tau_xy  =mu_previous.*(dudy_fwd+dvdx_central); % overwrite Tau_xy to change bias
    Tau_yy = 2*mu_previous.*(dvdy_fwd -1/3*(dudx_central+dvdy_fwd));
    qy = -k_previous.*dTdy_fwd;

    % 9: Create Fbar and take Derivative OF FBAR(backwards in y)
    F_bar(1,:,:) = U_bar(3,:,:);
    F_bar(2,:,:) = rho_previous.*u_previous.*v_previous - Tau_xy;
    F_bar(3,:,:) = rho_previous.*v_previous.^2+P_previous-Tau_yy;
    F_bar(4,:,:) = (Et_previous+P_previous).*v_previous-u_previous.*Tau_xy-v_previous.*Tau_yy+qy;

    for j = 1:4
       dFbardy_bwd(j,:,:) = ddy_bwd(squeeze(F_bar(j,:,:)),dy);
    end

    % 10: Calculate U_next using Ubar Fbar and Ebar (using vector column NS
    % eqn)
    U_next = 1/2*(U_previous+U_bar-dt*dEbardx_bwd-dt*dFbardy_bwd);

    % 11 - re - enforce BCs on all primitives before plotting

    % 11.1 - Convert U_next to prim vars (new Variables this time for some
    % reason?) (change this later)

    [rho,u,v,T,P,e,Et] = cons2prim(U_next,R,cv);
    a = real(sqrt(gamma*R*T)); %  recalculate local speed of sound at each point (at corrector level)
    mu = sutherland(T); % recalculate sutherland's law viscosity
    k = cp/Pr*mu; % recalculate conductivity with new primitives

    % 11.2.1 BCs - WALL & leading edge  (First may be overwritten by leading edge BC)
    T(:,1) = T_0;
    P(:,1) = extrapolate2(P,'bottom'); % 2nd order extrap. @ wall
    P(1,1)= P_0;
    u(:,1) = 0; % apply no-slip
    u(1,1) = u_0;
    v(:,1) = 0; % solid wall condition
    rho(:,1) = P(:,1)./(R*T(:,1)); % set rho BCs using ideal gas law
    rho(1,1) = rho_0 ; % This BC is not explicitly stated, but we must assume a BC for rho if P and T are fxed and P = rho * R *T
    e(:,1) = cv*T(:,1); % set sp. internal energy bcs using calorically perfect assumption
    Et(:,1) = rho(:,1).*(e(:,1)+1/2*u(:,1).^2);
    mu(:,1) = sutherland(T(:,1)); % recalculate viscosity with Temp BC in mind. 
    k(:,1) = cp/Pr*mu(:,1); % recalculae therm. conductivity with BC values of mu in mind

    % 11.2.2 BCs - INLET Boundary (Will not Overwrite Wall)
    T(1,2:end) = T_0;
    P(1,2:end) = P_0;
    u(1,2:end) = u_0;
    v(1,2:end) = 0 ;
    rho(1,2:end) = P_0/(R*T_0);
    e(1,2:end) = cv*T_0;
    Et(1,2:end) = rho(1,2:end).*(e(1,2:end)+1/2*(u(1,2:end).^2+v(1,2:end).^2));
    mu(1,2:end) = sutherland(T(1,2:end));
    k(1,2:end) = cp/Pr*mu(1,2:end);

    %11.2.3 - Far Field BCs (same condition as inlet, but on top boundary)
    %(overwrites outlet BC at corner)

    T(2:end,end) = T_0;
    P(2:end,end) = P_0;
    u(2:end,end) = u_0;
    v(2:end,end) = 0;
    rho(2:end,end) = P_0/(R*T_0);
    e(2:end,end) = cv*T_0;
    Et(2:end,end) = rho(2:end,end).*(e(2:end,end)+1/2*(u(2:end,end).^2+v(2:end,end).^2));
    mu(2:end,end) = sutherland(T(2:end,end));
    k(2:end,end) = cp/Pr*mu(2:end,end);

    % 11.2.4 Outlet BCS (as specified by prompt) will not overwrite far
    % field or wall conditions!
    
    T(end,2:(end-1)) = extrapolate2(T,'right');
    P(end,2:(end-1)) = extrapolate2(P,'right');
    u(end,2:(end-1)) = extrapolate2(u,'right');
    v(end,2:(end-1)) = extrapolate2(v,'right');
    rho(end,2:(end-1)) = P(end,2:(end-1))./(R*T(end,2:(end-1)));
    e(end,2:(end-1)) = cv*T(end,2:(end-1));
    Et(end,2:(end-1)) = rho(end,2:(end-1));
    mu(end,2:(end-1)) = sutherland(T(end,2:(end-1)));
    k(end,2:(end-1)) = cp/Pr*mu(end,2:(end-1));
    
    % 11.3 - reintroduce primitives with BCs updated into U_next

    U_next = prim2cons(rho,u,v,T,cv);

    % 12 Update this iterations primitives as "last iteration" in
    % preparation for the next loop
    u_previous = u;
    v_previous = v;
    e_previous = e;
    Et_previous = Et;
    mu_previous = mu;
    k_previous = k; 
    a_previous = a;
    rho_previous = rho;
    T_previous = T;
    P_previous =P;

    U_previous = U_next;


    %% Plotting and visualisation 
    
    if (mod(n,15)==0)

    subplot(2,3,1)
    pcolor(xx,yy,real(rho_previous))
    shading interp
    axis equal tight 
    c = colorbar;
    ylabel(c, 'Kg/m^3')
    xlabel('x');
    ylabel('y');
    title('\rho');

    subplot(2,3,2)
     pcolor(xx,yy,real(u_previous))
    shading interp
    axis equal tight 
    c = colorbar;
    ylabel(c, 'm/s')
    xlabel('x');
    ylabel('y');
    title('u');

    subplot(2,3,3)
    pcolor(xx,yy,real(v_previous))
    shading interp
    axis equal tight 
    c = colorbar;
    ylabel(c, 'm/s')
    xlabel('x');
    ylabel('y');
    title('v');

    subplot(2,3,4)
    pcolor(xx,yy,real(e_previous))
    shading interp
    axis equal tight 
    c = colorbar;
    ylabel(c,'J/kg')   
    xlabel('x');
    ylabel('y');
    title('e');

    subplot(2,3,5)
    pcolor(xx,yy,real(P_previous))
    shading interp
    axis equal tight 
    c = colorbar;
    ylabel(c, 'Pa')
    xlabel('x');
    ylabel('y');
    title('P');
    
    subplot(2,3,6)
    pcolor(xx,yy,real(T_previous))
    shading interp
    axis equal tight 
    c = colorbar;
    ylabel(c, 'K')
    xlabel('x');
    ylabel('y');
    title('T');

    sgtitle(sprintf('iteration Number %i',n));
    drawnow
   %% gif
    end

end



%% Functions in Use

function [mu] = sutherland(T)
mu_0 = 1.735*10^(-5); % Ns/m^2
T_0 = 288.15; 
S_1 = 110.4;

mu = mu_0*(T./T_0).^(3/2).*(T_0+S_1)./(T+S_1);
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


function boundaryVec = extrapolate2(F,boundary)
%% Extrapolate2 takes the 2nd order extrapolation of a provided field F (mxn) and boundary 'bottom' or 'right'
% if the right boundary is chosen, extrapolate2 returns a column vector of
% size mx1 which is the new boundary values of F, calculated via 2nd order
% extrapolation.  if the bottom boundary is chosen, a 1xn vector is
% returned. NOTE: THE FULL MATRIX F MUST BE PASSED TO extrapolate2, including boundary values which will then be overwritten!
arguments
    F
    boundary
end

[nx,ny] = size(F);

switch boundary 

    case 'right' % if user wants BCs for the rightmost side of F, apply 2nd order extrapolation there
    boundaryVec = zeros(1,(ny-2));
    for i=(1:ny-2) % NOTE: DOESNT COMPUTE EXTRAP @ top and bottom! Saves these vars for other bc calulations
        boundaryVec(1,i) = 2*F(end-1,i+1)-F(end-2,i+1); % 2nd order extrap. formula
    end

    case 'bottom' % if they want bottom, do the same there 
        boundaryVec = zeros(nx,1);
        for j=1:nx
            boundaryVec(j,1) = 2*F(j,2)-F(j,3);
        end
end

end
