format long
clear all 
close all 
clc


%% Define variables
Nb = 81; % Bin count

N_0 = 100; % simulated initial count of neutrons (ARBITRARY WILL WORK WITH ANY N0)
t_steps = 1000; % How many iterations of time steps should be calculated
q_star = .000735;%Neutron multiplcation factor for steady state peration (From HW3)
q = 3*q_star ; % as instructed, triple q for supercritical operation5555
a_plus = .873; % Reproduction rate inside control rod (slightly supercritical) .87282
a_minus= 0.872; %(subcritical alpha)
a_vec = [a_minus,a_plus];
Lbinloss = zeros(length(t_steps),1); % Tracks how many neutrons exit the left hand side of the reactor (stores a value for each time step)


%% Loop calculates and plots neutron counts in each bin as time progresses

%% PLEASE NOTE: BECAUSE OF LIVE UPDATING FIGURES,THIS LOOP TAKES ~2 MINS TO EXECUTE. 4 SUBPLOTS WILL GENERATE AFTER SUCCESSFUL COMPILATION
for m = 1:2 % Loop between both alpha values
    a = a_vec(m);
    % define in-loop arrays
    max_cr = zeros(t_steps,1); % Tracks the maximum neutron count in the critical reactor 
    max_cr(1) = N_0;
    bincount = zeros((t_steps),Nb+2); % creates a row for each time step, a column for each bin, 2 extra boundary bins which always have 0 neutrons
    bincount(1,2:end-1) = N_0; % sets initial values
    sum_n = zeros(t_steps,1); % records how many neutrons are present
    sum_n(1) = sum(bincount(1,:)); % Takes the first sum
  
    for i = 2:(t_steps) % loop through every time step
      Lbinloss(i-1,1)= 1/2*bincount(i-1,2);
    
          for n = 2:Nb+1 % loop through each bin (exluding the left and right boundary ones)
             bincount(i,n) = 1/2*bincount(i-1,n-1)+1/2*bincount(i-1,n+1);  % distribute half of the left bin into it, then half of the right bin
             
             if (i>1500)&&(n==42) % insert rod at 1500 Time steps, only in bin 42 (middle)eeeee
                 bincount(i,n) =bincount(i,n)*(a);
             else
                bincount(i,n) = bincount(i,n)*(1+q); % Multiply each bin by factor (1+q)
             end
    
           end
    
      sum_n(i,1) = sum(bincount(i,:)); % I just used this as a tool to track stability of the neutron count, not necessary for answers 1-3
      max_cr(i) = max(bincount(i,:));
    end

    if m==1
        max_alpha.subcritical = max_cr;
    else
        max_alpha.supercritical = max_cr;
    end
end

%% Plotting Neutron Count over time (to track stability)
% figure(1)
% 
% plot(1:t_steps,sum_n,'--b', Linewidth =2) % tracks overall neutron count with time
% xlabel('time step');
% ylabel('Total Neutrons present in reactor')
% title('Total Neutron Count in Reactor Over Time');
% xline(1500,'--r')
% legend('Total Neutron Count','Time of Control Rod Insertion')



  %% Plotting steady state neutron distribution

  subplot(2,1,1)
  plot(1:Nb,bincount(i,2:end-1),'.'); % Plots normalized neutron count
  title(' Time Evolution of Total Neutron Count Per Box')
  xlabel('Box Number')
  ylabel('Total Neutron Count');
  xlim([0,Nb+1]);

  subplot(2,1,2)
  plot(1:Nb,bincount(i,2:end-1)/max(bincount(i,2:end-1)),'.'); % Plots normalized neutron count
  title(' Time Evolution of Normalized Neutron Count Per Box (N/Nmax)')
  xlabel('Box Number')
  ylabel('N/Nmax');
  xlim([0,Nb+1])

 %% Question 1 plot: Tracking Nmax


figure(2)
hold on
time = 1:t_steps;
plot(time,max_alpha.supercritical,'.r') ;
plot(time,max_alpha.subcritical,'.b');
xline(1500,'--k');

xlabel('Timesteps')
ylabel('Nmax');
title('Maximum Neutron Count in the reactor for N_0 = 100');
legend('Neutron Count α_+ = 0.873', 'Neutron Count α_- = 0.872','Control Rod Insertion Time',Location='best')
hold off
