%% Q1


mu_e = 3.986*10^5;
v = 6.43; % km/s
r = 15000; % km
e = 0.825; 

        E = v^2/2-mu_e/r;
        a = -mu_e/(2*E);
        P = a*(1-e^2);
        h = sqrt(P*mu_e);

        phi = acosd(h/(r*v));






