%% Part 3: McCormack Advection 

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
T_bar = zeros(size(T));
T_halftime = zeros(size(T)); % stores the mean value of T n+1 and T n in each loop

figure(1)

%gif('Advection_McCormack.gif');

while current_time<2


    dTdx_fwd = ddx_fwd(T_new,dx,'p');
    dTdy_fwd = ddy_fwd(T_new,dy,'p');

    % first calculate the matrix of predicted T values T_bar for the next
    % timestep
    T_bar = T_new+dt_safe*(-1*c(1)*dTdx_fwd-c(2)*dTdy_fwd);

    % now calculate T new from Tbar
    dTbar_dx_bwd = ddx_bwd(T_bar,dx,'p');
    dTbar_dy_bwd = ddy_bwd(T_bar,dy,'p');
    % baked into the following calculation is the half timsetep
    % approximation of T n+1/2

    %% OLD CALLCULATION
    %T_new = T_new+dt_safe*(-1*c(1)*dTbar_dx_bwd-1*c(2)*dTbar_dy_bwd);
    %% NEW CALCULATION
    % note that the conversion from T_new was wrong previously because i
    % approximated T n+1/2 incorrectly.
    
    T_new = 1/2*((T_new+T_bar)+dt_safe*(-1*c(1)*dTbar_dx_bwd-1*c(2)*dTbar_dy_bwd));

    %plotting
    pcolor(xx,yy,T_new);
    colorbar
    clim([0,1]);
    axis equal tight
    shading interp
    xlabel('x')
    ylabel('y')
    title('Temperature - McCormack Method');
    
    subtitle_text = sprintf('Time =  %7.6f seconds after start',current_time);
    subtitle(subtitle_text);
    
   % gif
    current_time = current_time+dt_safe;
    drawnow

end