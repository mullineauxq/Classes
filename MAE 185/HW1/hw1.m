clear;
clc;
load('cylinder_Re100.mat');

%% Part 1 - Animated progression of the flow field 
t_steps  = length(u);

%commented gif command
% gif('dev_flow.gif')

for i = 1:t_steps
    subplot(2,1,1)
    pcolor(x,y,squeeze(u(i,:,:)));
    rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]);
    % add color scale limiting
    caxis([-0.5 2])
    title('u')
    xlabel('x')
    ylabel('y')
    axis equal tight
    shading interp

    subplot(2,1,2);
    pcolor(x,y,squeeze(v(i,:,:)));
    rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]);
    caxis([-1 1]);
    title('v')
    xlabel('x')
    ylabel('y')
    axis equal tight
    shading interp
    drawnow

   %gif
end
%% Part 2 - Mean flow field 

% retain, in a new MTRX, the last 151 values of u and v, removing the first
% 150, corresponding to a transient period in the flow field.
u_bar  = u(151:301,:,:); 
v_bar  = v(151:301,:,:);

%take the mean of each element in u and v (individually) across all 151
%'frames', then squeeze into a 2D array
u_bar_mean = squeeze(mean(u_bar));
v_bar_mean = squeeze(mean(v_bar));

% plot ubar mean, vbar mean 

subplot(2,1,1)
hold on

pcolor(x,y,u_bar_mean);
% fixed title
title('u mean');
xlabel('x');
ylabel('y');
shading interp
rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]);
axis equal tight
hold off
% added colorbar to emphasize scale
colorbar

subplot(2,1,2)
hold on
pcolor(x,y,v_bar_mean);
shading interp
title('v mean');
xlabel('x');
ylabel('y');
rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]);
axis equal tight
hold off
colorbar


%% Part 3 - Mean Flow Field and streamlines

% first, re-use the plot code from part 2 to generate the color map of the
% flow field "u". 

pcolor(x,y,u_bar_mean);
title('mean u');
xlabel('x');
ylabel('y');
shading interp
rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]);
axis equal tight
colorbar
hold off 

% I have no idea why, but this doesnt work unless I transpose everything...
    % transposing switches from ngrid to meshgrid!
% add more streamlines to the plot
streamline(x',y',u_bar_mean',v_bar_mean',-5*ones(14,1),[-4 -3 -2 -1.5 -1 -0.5 -0.01 0.01 0.5 1 1.5 2 3 4])
% add cylinder
rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]); hold off

%% Part 4 - Velocity Fluctuations 

% generall approach- loop through each u and v and subtract out the average
% developed ubar and vbar values 

gif('dev_flow_fluctuations.gif')
for i = 1:t_steps
    subplot(2,1,1)
    % express fluctuating component of u as u value minus ubar mean
    u_fluctuate_i = squeeze(u(i,:,:))- u_bar_mean;
    pcolor(x,y,u_fluctuate_i);
    rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]);
    % add plot accessories
    title('u fluctuation')
    xlabel('x')
    ylabel('y')
    axis equal tight
    colorbar
    shading interp

    subplot(2,1,2);
    v_fluctuate_i = squeeze(v(i,:,:))-v_bar_mean;
    pcolor(x,y,v_fluctuate_i);
    rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]);
  
    title('v fluctuation')
    xlabel('x')
    ylabel('y')
    axis equal tight
    shading interp
    colorbar
    drawnow

   gif
end
