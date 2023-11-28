%% Define the system transfer function with time delay
numerator = [.1];
denominator = [ 1, .1]; 
time_delay = 6; %  time delay in seconds

g = tf(numerator, denominator);
g_delayed = tf(numerator,denominator, 'InputDelay', time_delay);
g_delayed_rational = pade(g_delayed,2);


%% Create a Bode plot
figure(1);
hold on
bode(g);
bode(g_delayed)
grid on;
title('Bode Plot with Time Delay');
legend("No delay","delayed 6s"); 
hold off
%% 

%%d = tf([ 4.94 2 0.20249],[1 0]); %% PID controller found using Zeigler
%%Nichols (OLD using 2pi factor)

d = tf([.786 2 1.274],[1 0]);
L = g*d ; %% Open loop transfer func
L_delayed = g_delayed*d; 

Ts = L_delayed/(1+L_delayed); %% configure closed loop tf T(s)
[numerT,denomT] = tfdata(Ts);

figure(2)
hold on
bode(L)
bode(L_delayed)
title('Bode Plot of Forward Loop using Ziegler Nichols PID');
hold off
legend("No delay","delayed 6s"); 

%% 4a Plotting System response to sequential step inputs

t = [0:1:18000];
u = zeros(1,18001);
u(1)=20;
for n=2:3600
u(n) = 35;
end
for n=3601:14401
u(n) = 45;
end
for n=14401:18001
u(n) = 20;
end

figure(3)


y = lsim(Ts,u,t);
ylabel('Temperature C');

[Gm,Pm,Wcg,Wcp] = margin(L_delayed);



