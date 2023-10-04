lam_p = 0.0001 ;% simulated value of lambda P (parent isotope decay constant (1/s)
lam_d= 1000*lam_p; % must be larger significantly than parents decay constant 

Np0 = 10000;
Nd = @(t) lam_p*Np0/lam_d*(1+exp(-lam_d*t));

figure(1)
hold on
fplot(Nd,'-b',LineWidth=2);

yline(lam_p*Np0/lam_d,'--r',LineWidth=2);

xlim([0,100])

ylim([0 20]);

xlabel('Time');

ylabel('Atoms of Daughter Isotope');

title('Secular equilibrium of Daughter isotope');

legend('Daughter Isotope count','Equilibrium value = lamP*Np0/lamD')