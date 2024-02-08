%% Code for reading the ECP codes, geenerating the plot and getting the figures of interest

ECP_vals  = readecp('OL_3_configurationB.txt');

time = ECP_vals(:,2);

encoder2 = ECP_vals(:,5);

figure(1)

plot(time,encoder2)

title('Encoder 2 Step response');

xlabel('Time (s)');

ylabel('Position (counts)');


