N = [1 1];
D= [1 10 1];
Ts = 0.1;
f=  0;

% Takes coefficeint vector of numerator N and denominator D as a symbolic matrix with 2 rows. 1st row vector of the
% numerator coefficients, second of the denominator coefficients.
% Take the symbolic s-domain tf represented by these values and create a
% matched Z transform Gz using sample time Ts, and target frequency f (defaults to 0). Set
% the causal field to 1 if strictly causal Gz is required, 0 if not. 

h = Ts;
f = f;
zeros = roots(N);
m = length(zeros);
poles = roots(D);
n = length(poles);
v = n-m; % used to add infinite zeros in s to z domain 



Dz = 1; % initialize z transfer function
syms z

for i =1:length(zeros)
Dz = Dz*(z-exp(zeros(i)*h));
end

for j=1:length(poles)
Dz = Dz*1/(z-exp(poles(j)*h));
end

if v>2
Dz = Dz*(z+1)^(v-1);
end


sympref('FloatingPointOutput',true);
Dz_matlab = c2d(tf(N,D),0.1)


%% Gain matching 
syms s
K_s = limit(poly2sym(N,s)/poly2sym(D,s),s,i*(f));

if nnz(f==poles)||nnz(f==zeros) % if the frequency is a pole or zero of the S domain TF, match a close frequency that isnt this one
    f = f+0.001; % slightly adjust f
    K_s = limit(poly2sym(N,s)/poly2sym(D,s),s,i*(f)); % if gain is zero, gain match at a different frequency close to desired 
end

K_z = limit(Dz,z,exp(i*(f)));

Keff = K_s/K_z;

Dz = Dz*Keff
