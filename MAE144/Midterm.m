%% Define the system transfer function with time delay
numerator = [.1];
denominator = [ 1, .1]; 
time_delay = 6; %  time delay in seconds

g = tf(numerator, denominator);
g_delayed = tf(numerator,denominator, 'InputDelay', time_delay);
g_delayed_rational = pade(g_delayed,2);

Ts_delayed = 3.4*g_delayed_rational/(1+3.4*g_delayed_rational);

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

d = tf([616.4 1568.4 999],[785 0]); %% PID controller found using Zeigler Nichols
L = g*d ; %% Open loop transfer func
L_delayed = g_delayed_rational*d; 

Ts = L/(1+L); %% configure closed loop tf T(s)
[numerT,denomT] = tfdata(Ts);

T_delayed = tf(numerG,denomG,'InputDelay',6);

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


[NUM,DEN] = tfdata(T_delayed);
Num = NUM{:};
Den = DEN{:};

[A,B,C,D] = tf2ss(Num,Den);

T_delayed_ss = ss(A,B,C,D);
lsim(T_delayed,u,t);
ylabel('Temperature C');

[Gm,Pm,Wcg,Wcp] = margin(T_delayed);


