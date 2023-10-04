format long
clear all 
close all 
clc


%% Define variables
Nb = 81; % Bin count
N_0 = 100; % simulated initial count of neutrons (ARBITRARY WILL WORK WITH ANY N0)
t_steps = 1000; % How many iterations of time steps should be calculated

q = .000735; %Neutron multiplcation factor which I tuned to reach a steady state

Lbinloss = zeros(length(t_steps),1); % Tracks how many neutrons exit the left hand side of the reactor (stores a value for each time step)

%% Define arrays and structures


bincount = zeros((t_steps),Nb+2); % creates a row for each time step, a column for each bin, 2 extra boundary bins which always have 0 neutrons
bincount(1,2:end-1) = N_0; % sets initial values

sum_n = zeros(t_steps,1); % records how many neutrons are present
sum_n(1) = sum(bincount(1,:)); % Takes the first sum

max_cr = zeros(t_steps,1); % Tracks the maximum neutron count in the critical reactor 
max_cr(1) = N_0;

%% Loop calculates and plots neutron counts in each bin as time progresses

%% PLEASE NOTE: BECAUSE OF LIVE UPDATING FIGURES,THIS LOOP TAKES ~2 MINS TO EXECUTE. 4 SUBPLOTS WILL GENERATE AFTER SUCCESSFUL COMPILATION

for i = 2:(t_steps) % loop through every time step
  Lbinloss(i-1,1)= 1/2*bincount(i-1,2);

      for n = 2:Nb+1 % loop through each bin (exluding the left and right boundary ones)
         bincount(i,n) = 1/2*bincount(i-1,n-1)+1/2*bincount(i-1,n+1);  % distribute half of the left bin into it, then half of the right bin
         bincount(i,n) = bincount(i,n)*(1+q); % Multiply each bin by factor (1+q)
      end 

  sum_n(i,1) = sum(bincount(i,:)); % I just used this as a tool to track stability of the neutron count, not necessary for answers 1-3
  max_cr(i) = max(bincount(i,:));
  %% Plotting Question 1 
  subplot(3,1,1)
  plot(1:Nb,bincount(i,2:end-1),'.'); % Plots normalized neutron count
  title(' Time Evolution of Total Neutron Count Per Box')
  xlabel('Box Number')
  ylabel('Total Neutron Count');
  xlim([0,Nb+1]);

  subplot(3,1,2)
  plot(1:Nb,bincount(i,2:end-1)/max(bincount(i,2:end-1)),'.'); % Plots normalized neutron count
  title(' Time Evolution of Normalized Neutron Count Per Box (N/Nmax)')
  xlabel('Box Number')
  ylabel('N/Nmax');
  xlim([0,Nb+1])
  %drawnow
end



%% Plotting Other Data (Helps track Stability of reactor)
% figure(2)
% 
% plot(1:t_steps,sum_n,'--b', Linewidth =2) % tracks overall neutron count with time
% xlabel('time step');
% ylabel('Total Neutrons present in reactor')
% title('Total Neutron Count in Reactor Over Time');


%% Question 2- Normalize left side loss with Jdiff

Nmax = max(bincount(end,:)); % Finds the highest single bin neutron count across all bins
jdiff_dt = Nmax*(pi)/(2*Nb); % This is the analytical formula for Jdiff*dt I found.

subplot(3,1,3)
hold on
plot(1:t_steps-1,Lbinloss/jdiff_dt, '-r', linewidth=2);
title('Normalized Neutron Loss Through the Left Boundary per Timestep');
xlabel('Time Step');
ylabel('Nloss/Jdiff*dt')

hold off

%% Question 3 


q_star = 4*q ;

bincount_sc = zeros(t_steps,Nb+2);% Bincount of neutrons in the now supercritical reactor
max_sc = zeros(length(t_steps),1); % tracks the maximum number of Neutrons in any bin

bincount_sc(1,2:(end-1))=N_0; % set first bin values 
max_sc(1)= max(bincount_sc(1,2:(end-1))); % record maximum before first iteration

for i = 2:(t_steps) % loop through every time step

      for n = 2:Nb+1 % loop through each bin (exluding the left and right boundary ones)
         bincount_sc(i,n) = 1/2*bincount_sc(i-1,n-1)+1/2*bincount_sc(i-1,n+1);
         % distribute half of the left bin into it, then half of the right bin
         bincount_sc(i,n) = bincount_sc(i,n)*(1+q_star); % bin multiplier
      end 

max_sc(i,1) = max(bincount_sc(i,:)); % record maximum for this time step

end

time = 1:t_steps;% created only for plotting purposes

figure(3)

semilogy(time,max_sc); % Data-driven result curve 

hold on

N_o = max_sc(end)/exp((q_star-pi^2/(2*Nb^2))*t_steps); % backwards - calculate n_0 (Yutong's method)

n_t = @(t) N_o*exp((q_star-pi^2/(2*Nb^2))*t); % implicit function for n(t), time dependent solution for neutron distribution


semilogy(time,n_t(time)); %plot analytical result

legend('Maximun Neutron Count', 'Analytical Solution')
xlabel('Time Steps');
ylabel('Log Maximum N count');
title('Maximum Neutron Count over Time Compared to Analytical');

%% Track Nmax over time to prove stability 

figure(4)

plot(time, max_cr,'.r')
xlabel('Timestep')
ylabel('Nmax');
title('Tracking Maximum Neutron Count');




