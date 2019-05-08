% author : Ashwin de Silva - 150103P

%% Schrodinger Wave Equation Solutions in 1D for Pre-Defined Potential Functions - NMM

clear all;
close all;

%% constants
h_bar = 1.055e-34;
m_e = 9.11e-31;
e = 1.602e-19;

%% inputs 
Lz = 10e-10;    % width of the potential well
N = 1000;       % Number of z values
Vo = 30;
Wz = Lz/4;

%% setting up z coordinates
z_min = 0;
z_max = Lz;
z = linspace(z_min,z_max,N);
dz = z(2)-z(1);

%% define the potential

% f = @(x) Vo/(Lz)*x;
% f = @(x) Vo/Lz^2*x.^2;
% f = @(x) Vo/(Wz-Lz)^2*(x + Lz).^2; 
% f = @(x) Vo/(Lz-Wz)*x + Vo + Vo/(Lz-Wz)*Wz;
% f = @(x) Vo;

% V = potential_func(z,Wz,f,1)'*e;

% V = stepped_well(z,Wz,Vo)'*e;
% V = double_well(z,Wz,Vo)'*e;
% V = assymetric_well(z,Wz,Vo)'*e; 
V = sloping_well(z,Wz,Vo)'*e;
% V = morse_potential(z,Wz,Vo)'*e;

%% define the matrices

i = ones(N,1); 
Lap = spdiags([i -2*i i],[-1 0 1],N,N)/dz^2;

%% create the Hamiltonian

H = -1/2*(h_bar^2/m_e)*Lap + spdiags(V,0,N,N);

%% find the eigenfunctions and eigenvalues of the H
nmodes = 3; options.disp = 0;
[Vec,Eig] = eigs(H,nmodes,'sa',options); % find eigs
[Eig,ind] = sort(diag(Eig));
Vec = Vec(:,ind);

%% Display results

disp('Allowed Energies (eV) : ');
disp(Eig/e);

numerov_plotter(Vec,V/e,z);
energy_plotter(Eig/e,z,V/e)





