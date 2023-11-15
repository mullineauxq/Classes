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