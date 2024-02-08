%% Problem1
clc 
clear

load('gotritons','T','xx','yy')
delta_x = mean(diff(xx(:,1))); 
delta_y = mean(diff(yy(1,:)));
alpha = 2; %a given value
factor_of_safety = 1.2 ;

%% calculate the minimum timestep using stability condition
dt_stable = 1/(alpha*4)*(delta_x^2*delta_y^2)/(delta_x^2+delta_y^2); % analytical stability requirement
dt_safe = dt_stable/factor_of_safety; % calculate timestep



T_new = T; % tracks the value of T over time
delta_T = zeros(400,200); % tracks the chages in T in each time step 
t_current= 0;
figure(1)
%gif('temperature_diffusion.gif') % this line saves a gif only with a
%plugin 

while t_current <0.001
% calculate change in temp assumming forward first diff in time, periodic
% 2nd diff in space
delta_T = dt_safe*alpha*(d2dy2(T_new,delta_y,'p')+d2dx2(T_new,delta_x,'p'));

% add dT to the temp from the last iteration
T_new = T_new+delta_T;

pcolor(xx,yy,T_new);
colorbar
clim([0,1]);
axis equal tight
shading interp
xlabel('x')
ylabel('y')
title('Temperature over test area');

subtitle_text = sprintf('Time =  %7.6f seconds after start',t_current);
subtitle(subtitle_text);

t_current = t_current+dt_safe;
%gif
drawnow

end



