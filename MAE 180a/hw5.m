format long
r1 = [1.81065659, 1.06066883, .31065150];
r2 = [1.41422511,0,1.414202];
r3 = [1.35353995, 1.41422511, -0.64644950];

D = cross(r1,r2)+cross(r2,r3)+cross(r3,r1);
d = norm(D);

N = norm(r3)*cross(r1,r2)+norm(r1)*cross(r2,r3)+norm(r2)*cross(r3,r1);
n=norm(N);

S = (norm(r2)-norm(r3))*r1+(norm(r3)-norm(r1))*r2+(norm(r1)-norm(r2))*r3;
s = norm(S);

L = sqrt(1/(d*n));

v1 = L/norm(r1)*cross(D,r1)+L*S;

h = cross(r1,v1);

e =cross(v1,h)-r1/norm(r1);

E = norm(v1)^2/2-1/2.12;

a = -1/(2*E);

P = a*(1-norm(e)^2);

i= acosd(dot(h,[0,0,1])/norm(h));

node = cross([0 0 1], h);

omega = acosd(dot(node,e)/(norm(node)*norm(e)));
OMEGA = acosd(dot(node,[1 0 0])/(norm(node)));

nu_o1 = acosd(dot(e,r1)/(norm(e)*norm(r1)));


%% question 2
el = 40;
az = 65;
rho = .8*[-cosd(el)*cosd(az);cosd(el)*sind(az);sind(el)];

vr = [-.07*cosd(40)*cosd(65)+.8*2.5*sind(40)*cosd(65)+0.8*cosd(40)*1.8*sind(65);
    0.07*cosd(el)*sind(az)-0.8*2.5*sind(el)*sind(az)+0.8*cosd(el)*1.8*cosd(az);
    0.07*sind(el)+0.8*2.5*cosd(el)];