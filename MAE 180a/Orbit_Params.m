function Params =  Orbit_Params(r,v, mu)
% This function takes in a symbolic or numeric vector v and a symbolic or
% numeric vecotr r of the same size. r and v, describing the orbit radius
% and velocity, can be used to calculate all of the orbital parameters h (angular momentum vector),
% e(eccentricity vector), i (inclination angle btw z axis and h), LAN
% ( longitude of the ascending node  (from X axis)), omega (ascending node to periapsis
% angle), and nu_0 (true anomaly). mu defaults to earth's value of mu in
% km^3/s^2. set it otherwise if using a different orbiting body or units.

arguments
    r; % no default value 
    v;
    mu = 3.986*10^5; %km^3/s^2 , the default auto sets to earth in SI units
end

h = cross(r, v);
Params.h = h;

e = 1/mu*((norm(v)^2-mu/norm(r))*r-(r*v')*v);
Params.e = e;

k = [0 0 1];
i = acos(h*k'/norm(h));
Params.i = i ;

 

n = cross(k,h); % define the ascending node direction
i_hat = [1 0 0];
OMEGA = acos(n*i_hat'/(norm(n)));

if n(2) < 0 
    OMEGA = 2*pi-OMEGA; % condition for ascending node in the 3rd or 4th quadrant of the equatorial plane
end
Params.LAN = OMEGA ; 

omega = acos(n*e'/(norm(n)*norm(e)));

if e(3)<0
    omega = 2*pi-omega; % CONDITION FOR PERIAPSIS DIRECTION BELOW EQUATORIAL PLANE (-Z)
end

Params.omega = omega;


nu = acos(e*r'/(norm(e)*norm(r)));

if r*v' < 0 
nu = 2*pi-nu;
end

Params.anomaly = nu;

end