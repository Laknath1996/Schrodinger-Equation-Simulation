% author : Ashwin de Silva - 150103P

%% Schrodinger Wave Equation Solutions in 1D for Pre-Defined Potential Functions - FDM

clear all;
close all;

%% constants
h_bar = 1.055e-34;
m_e = 9.11e-31;
e = 1.602e-19;

%% inputs 
Lz = 8e-10;                     % width of the well (m)
N = 1000;                       % resolution of the z axis
E1 = 0;                         % lower bound of the energing shooting (eV)
E2 = 10;                        % upper bound of the energing shooting (eV)
steps = 100000;                 % steps resolution in energy shooting

%% get the energy steps
E = E1:1/steps:E2;              
dE = E(2) - E(1);

%% setting up z coordinates
z_min = -Lz;
z_max = Lz;
z = linspace(z_min,z_max,N);
dz = z(2)-z(1);

%% potential function

Vo = 10;                                        % highest potential
Wz = Lz/2;                                      % well boundary

% define the function

% f = @(x) Vo/(Lz)*x;                           % Infnite Well with Field
% f = @(x) Vo/Lz^2*x.^2;                        % Harmonic Oscillator
% f = @(x) Vo/(Wz-Lz)^2*(x + Lz).^2;            % Quadratic Potential Well
% f = @(x) Vo/(Lz-Wz)*x + Vo + Vo/(Lz-Wz)*Wz;   % Linear Potential Well
% f = @(x) Vo;                                  % Square Potential Well
% f = @(x) Vo/Wz^2*x.^2;                        % Truncated Parabolic Well
% f = @(x) abs(Vo/Wz*x);                        % Truncated |x| well


V = potential_func(z,Wz,f,2);                   % for potential wells : 1, for non potentials wells : 0
% V = zeros(1,length(z));

%%  Find the eigen-energies by Finite Difference Method

nozc = zeros(1,length(E));                      % non zero count
tv = zeros(1,length(E));                        % terminal values

for m = 1:length(E)                             % energy shooting
    psi = wav_func(E(m),V,dz,N);
    zc = zero_cross_count(psi);
    nozc(m) = zc;
    [zc,id] = zero_cross_count(psi);
    if range(id(1:end-1)) - range(id(2:end)) < 2
        nozc(m) = zc;
    else
        nozc(m) = -1;
    end
    tv(m) = psi(N);
end
E(nozc == -1) = [];
tv(nozc == -1) = [];
nozc(nozc == -1) = [];
tv = abs(tv);

% Determine the number of allowed eigen states in the provided energy range
es = unique(nozc);

% check for the higest eigen-state in the provided energy range 
ee = zeros(1,length(es));    
    
%% Find the eigen-energies for each eigenstate
for i = 1:length(es)
    [mtv,idx] = min(tv(nozc == es(i)));
    E_temp = E(nozc == es(i));
    ee(i) = E_temp(idx);
end

%% Display the results

disp(sprintf('Number of Eigen States between %d eV and %d eV : %d',E1,E2,length(es)));
disp('Allowed Energies (eV) : ');
disp(ee);

plotter(ee,z,V,dz,N);
energy_plotter(ee,z,V);
