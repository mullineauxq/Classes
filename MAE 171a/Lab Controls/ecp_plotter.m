%% Plot encoder 2 position, commaned positionover time


Data = readecp('2DOF_CL_exp2.txt'); %Change this filename to your file!

time = Data(:,2);

commanded_position = Data(:,3);

encoder2_position = Data(:,5);

figure(1)

hold on
plot(time, commanded_position,LineWidth=2);
plot(time, encoder2_position,LineWidth = 2);
legend('Commanded Position','Encoder 2 Position');

hold off


