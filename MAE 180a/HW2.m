%% Q1


mu_e = 3.986*10^5; %km^3/s^2
re = 6378.137; %EARH RADIUS KM
v = 6.43; % km/s
r = 15000; % km
e = 0.825; 


        E = v^2/2-mu_e/r;
        a = -mu_e/(2*E);
        P = a*(1-e^2);
        h = sqrt(P*mu_e);

        phi = acosd(h/(r*v));

 % part b

r2=110+re;
v2 = sqrt(2*(E+mu_e/r2));

    phi = acosd(h/(r2*v2));

    rp = a*(1-e);



%% problem 2

mu_s = 3.793119*(10^7);
rs = 60268; % saturn radius
Es = -mu_s/(2*rs);

vs = sqrt(2*(Es+mu_s/rs));
TU= rs/vs;
DU=rs;

r = [0 0 1.4];
v=[sqrt(3) sqrt(3) 0];
h = cross(r,v);
phi = acosd(sqrt(h*h')/(sqrt(r*r')*sqrt(v*v')));

E = norm(v)^2/2-1/1.4;

e = sqrt(1+2*E*norm(h)^2);

%% Problem 3

rp = 140 +1737.5 ; % radius of satellite above moon km
mu_m = 4.902801*10^3; %mu moon 

dv = sqrt(mu_m/rp)*(sqrt(2)-1);

dve = sqrt(2*6.8450+2*mu_m/rp)-sqrt(mu_m/rp);


%% problem 4
mu_e = 3.986*10^5; %km^3/s^2
E = 3.916;
e = 1.207;
a = -mu_e/(2*E);

rp = a*(1-e);

Ec = -mu_e/(2*rp);

vc= sqrt(2*(Ec+mu_e/rp));

vp = sqrt(2*(3.916+mu_e/rp));

dv= vc-vp;




