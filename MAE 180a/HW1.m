format long
mu = 3.986*10^5;
r = [8228 389 6888];
v = [-0.7 6.6 -0.6];

h = cross(r,v);

norm(h);

E = v*v'/2 - mu/norm(r); % Energy of the Orbit

P =h*h'/mu ; % The orbit Parameter

e = cross(v,h)/mu - r/norm(r);

