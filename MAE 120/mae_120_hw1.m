
format long
%% Define Masses (amu) and constants
Mn = 1.0086649157; % Netutron mass (amu)
Mp = 1.00727646676; % Proton mass 
Me = 0.000548579911; % Electron
M_N14 = 14.0030740048;
M_O17 = 16.99913170; 
M_C12 = 12;
M_N13 = 13.00573861;
M_He4 = 4.00260325415;
M_U235 = 235.04392992;
MH = 1.00782503207;
M_H3 = 3.0160492777;
M_Li6= 6.015122795;


c = 299792458; % m/s
Na = 6.0221408*10^(23); % avogadros (amu/g)
K = 931.494 ; % MeV/amu

%% Atomic  Calculations

dm_He4 = M_He4-2*Mp-2*Mn-2*Me;
dm_U235 = 143*Mn + 92*MH - M_U235;
%binding energy in eV
re_He4 = dm_He4/(1000*Na)*c^2*6.241509*10^18; % E =mc^2 converted to kg
re_U235 = dm_U235/(1000*Na)*c^2*6.241509*10^18;

%% Problem 5

dm_rxn = M_Li6+Mn-M_He4-M_H3;

%% Problem 6 

dm_rx1 = M_He4 +M_N14 - M_O17-MH;
E1 = K*dm_rx1;
dm_rx2 = MH+M_C12-M_N13;
E2 = K*dm_rx2;
dm_rx3 = Mn +M_Li6 -M_H3 - M_He4;
E3 = K*dm_rx3;
