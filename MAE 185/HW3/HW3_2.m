%% Part 2: Euler Advection 

load('gotritons','T','xx','yy')

dx = mean(diff(xx(:,1))); 

dy = mean(diff(yy(1,:)));

% define terms  
c = [1;1];
safety_factor = 2;
dt_stable  = 1*(dy*dx/(dy+dx));
dt_safe = dt_stable/safety_factor; 

%% Compute time evolution with first order backwards diff until t=2
T_new = T; 
current_time = 0 ; % tracks the time with each iteration
dT = zeros(size(T));


figure(1)
%gif('advection_backwards diff.gif')

while current_time<2
    dTdy = ddy_bwd(T_new,dy,'p');
    dTdx = ddx_bwd(T_new,dx,'p');
    
    dT = dt_safe*(-c(1)*dTdx - c(2)*dTdy);

    T_new = T_new+dT;

    pcolor(xx,yy,T_new);
    colorbar
    clim([0,1]);
    axis equal tight
    shading interp
    xlabel('x')
    ylabel('y')
    title('Temperature');
    
    subtitle_text = sprintf('Time =  %7.6f seconds after start',current_time);
    subtitle(subtitle_text);
    
  %  gif
    current_time = current_time+dt_safe;
    drawnow
end

%% now attempt the same procedure but with a centered diff

current_time = 0;
T_new = T; 
figure(2)
gif('advection_center_diff.gif')
dT = zeros(size(T));

while current_time<0.25
    dTdy = ddy_central(T_new,dy,'p');
    dTdx = ddx_central(T_new,dx,'p');
    
    dT = dt_safe*(-c(1)*dTdx - c(2)*dTdy);

    T_new = T_new+dT;

    pcolor(xx,yy,T_new);
    colorbar
    clim([0,1]);
    axis equal tight
    shading interp
    xlabel('x')
    ylabel('y')
    title('Temperature (Center Diff Method)');
    
    subtitle_text = sprintf('Time =  %7.6f seconds after start',current_time);
    subtitle(subtitle_text);
    
    gif
    current_time = current_time+dt_safe;
    drawnow
end
