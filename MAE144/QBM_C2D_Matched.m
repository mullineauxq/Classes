function [Dz,Keff] = QBM_C2D_Matched(N,D,Ts,f,causal)
% Takes coefficeint vector of numerator N and denominator D as a symbolic matrix of polynomial coeffcieints of S with 2 rows. 1st row vector of the
% numerator coefficients, second of the denominator coefficients.
% Take the symbolic s-domain tf represented by these values and create a
% matched Z transform Gz using sample time Ts (defaults to 0.01), and target frequency f (defaults to 0). Set
% the causal field to 1 if strictly causal Gz is required, 0 if not.
% (defaults to causal). Returns both Dz (matched z transfer function) and
% Keff, the constaint gain factor to multiply against Dz to gain match it
% against D(s) at the chosen frequency f. 

arguments
    N
    D
    Ts = 0.010; % sets up default value to pass, allowing user to skip inputting this
    f = 0; % default gain-matching frequency
    causal= 1; % defaults to causal Z-domain transform 
end
    h = Ts;

    zeros = roots(N);
    m = length(zeros);
    poles = roots(D);
    n = length(poles); 
    v = n-m; % used to add infinite zeros in s to z domain


    Dz = 1; % initialize z transfer function
    syms z

    for i =1:length(zeros) % convert all zeros to Z domain by Z = e^(s*h)
        Dz = Dz*(z-exp(zeros(i)*h));
    end

    for j=1:length(poles) % same for all poles
        Dz = Dz*1/(z-exp(poles(j)*h));
    end

    if v>2  % if numerator order is v above denom, where v>2 , add in v-1 zeros at s=inf, z=-1
        Dz = Dz*(z+1)^(v-1);
    end

    if causal&& v>2 
        Dz = Dz/(z+1); % remove a zero at z=-1 if causality is required
    end

    sympref('FloatingPointOutput',true); % helps keep long rationals out of the returned answer
    %c2d(tf(N,D),0.1)


    %% Gain matching
    syms s
    K_s = limit(poly2sym(N,s)/poly2sym(D,s),s,i*f); % gain of CT TF at frquency f 

    if nnz(f==poles)||nnz(f==zeros) % if the frequency is a pole or zero of the S domain TF, match a close frequency that isnt this one
        f = f+0.001; % slightly adjust f
        K_s = limit(poly2sym(N,s)/poly2sym(D,s),s,i*(f)); % if gain is zero, gain match at a different frequency close to desired
    end

    K_z = limit(Dz,z,exp(i*(f))); % gain of the matched Z domain tf at frequency f

    Keff = K_s/K_z; % effective multiplication factor to multiply D(z) by in order to gain match it against D(s) at frequency f

    Dz = simplify(Dz); % may help remove further extraneous rationals

end