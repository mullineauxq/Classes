clear 
clc

syms s
%% Define G, T (Known Polynomials)
numer_g = (s+2)*(s-2)*(s+5)*(s-5); % Define G numerator
denom_g = (s+1)*(s-1)*(s+3)*(s-3)*(s+6)*(s-6); % define G denominator
N_g = RR_poly(sym2poly(numer_g)); % format polynomial (from numerator of G) into coefficient vector for RR_poly

D_g = RR_poly(sym2poly(denom_g));% coefficeint vector for G denominator converted to     RR_poly

denom_t = (s+1)^2*(s+3)^2*(s+6)^2; % Set up the polynomial described in prompt 2a
D_t = RR_poly(sym2poly(denom_t)); % Classify the feedback loop denominator as an RR_poly 


[D_d,N_d] = RR_diophantine(D_g,N_g,D_t); % Call the Diophantine function to give polynomial coefficients to D_g and N_g


D_d*D_g+N_g*N_d; %% Simple check that the resulting polynomial has the roots prescribed by the prompt



%% Problem  2b


% How large does k (additional repeat root poles @ s=-60) have to be to
% make D(s) proper?

m = length(N_d.poly); % Find the order of the numerator polynomial
n = length(D_d.poly); % order of the denominator polynomial
k=0;

while m>n % stops when m<=n (properness condition)
k=k+1;

denom_t = (s+1)^2*(s+3)^2*(s+6)^2*(s+20)^k; % Added a pole at s=-60 multiplicity k
D_t = RR_poly(sym2poly(denom_t)); % convert to RR_poly
[D_d,N_d] = RR_diophantine(D_g,N_g,D_t); % call Diophantine

m = length(N_d.poly); % Recalculate orders
n = length(D_d.poly); 
end

fprintf('It takes a repeat pole of multiplicty %1.0f in f(s) to make D(s) proper using the Diophantine code \n',k);





