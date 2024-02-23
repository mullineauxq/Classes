load('gotritons_uv.mat')


t = 0; % tracks loop time 
dt = 0; % tracks time step 
FoS = 1.2;% factor of safety for stable timestep choice
d = 0.0005; % diiffusion coeff as specified in problem statement

dx  = mean(diff(xx(:,1)));% average 'dx' spacing
dy = mean(diff(yy(1,:)));

u_last = u ; % tracks the values of u and v on the previous iteration
v_last = v ;

u_bar = zeros(size(u)); % tracks the predicted values for the next timestep 
v_bar = zeros(size(v));


u_next = zeros(size(u)); % tracks the true values for the next timestep 
v_next = zeros(size(v));

%each loop, user euler explicit step for u and v n+1, use the spatial
%derivatives as specified in the problem statement 
% Loop until 2s is reached (as specified in problem statement)

figure(1)
u_last(:,1) = 0 ; % x velocity 0 at top and bottom boundary
u_last(:,end) = 0;
%gif('x_velocity_burgers_eqn.gif');
while t<2 

%enforce boundary conditions(PERIODICITY HANDLED WITHIN derivative
%functions)

v_last(:,1) = 0;  % y velocity at top and bottom  = 0 
v_last(:,end) = 0; 
% end boundary condititons 

umax = max(max(abs(u_last)));% calculate maximum velocities from previous timestep
vmax = max(max(abs(v_last)));

% update timestep length and current time 
dt = t_safe_mccormack(umax,vmax,dx,dy,FoS);
t = t+dt;

% calculate forward and central differences 
% recall that the boundary is periodic only in x 

d2udx2 = d2dx2(u_last,dx,'p'); % central diffs for 2nd order, periodic only in x
d2udy2 = d2dy2(u_last,dy);

dudx_fwd = ddx_fwd(u_last,dx,'p'); % forward differences in space, pin xeriodic only 
dudy_fwd = ddy_fwd(u_last,dy);

%Using forward diffs in x y and t, calculate a predicted value for u in the
%next time step (v omitted due to zero state initial condition

u_bar = u_last + dt*(d*(d2udx2 +d2udy2)-u_last.*dudx_fwd-v_last.*dudy_fwd);

u_bar(:,1) = 0 ; % re-establish no-slip bcs
u_bar(:,end) = 0 ;

% with U bar calculated, call a euler explicity in t, centered 2nd diff in
% space and backwards 1st diff in space where applicable on u_bar.

d2ubardx2 = d2dx2(u_bar,dx,'p') ;% calculate those centered 2nd diffs on u_bar, periodic in x
d2ubardy2 = d2dy2(u_bar,dy);
dubardx_bwd = ddx_bwd(u_bar,dx,'p'); % periodic backwards 1st diff in x
dubardy_bwd = ddy_bwd(u_bar,dy);

% now calculate the true next u from the previous iterations' u_last, and
% the predictor u_bar.
u_next = 1/2*dt*(d*(d2ubardx2+d2ubardy2)-u_bar.*dubardx_bwd-v_last.*dubardy_bwd)+1/2*(u_bar+u_last);
u_next(:,1) = 0; % enfoce bcs before plotting
u_next(:,end) = 0 ;

pcolor(xx,yy,u_next)
c = colorbar;
shading interp
axis equal tight 
xlabel('x')
ylabel('y')
ylabel(c,'Velocity m/s')
title('Burgers Equation: X velocity',sprintf('t = %f seconds',t))
u_last = u_next;
drawnow
%gif
end
% end while loop


function dt = t_safe_mccormack(umax,vmax,dx,dy,FoS)
% calculates a safe timestep to use in mcCormack method
dt = dx*dy/(umax*dy+vmax*dx); % using Mccormacks dt formula
dt = dt/FoS; % add factor of safety
end