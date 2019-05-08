% author : Ashwin de Silva - 150103P

%% Schrodinger Wave Equation Solutions in 2D for Pre-Defined Potential Functions

clear all;
close all;

%% User Inputs
num_sol = 15;    % desired number of solutions
Vo = 15;       % desired potential barrier in electron volts  
Wx = 50e-10;    % x width of the rectangular well
Wy = 70e-10;    % y width of the rectangular well
Rx = 50e-10;    % x semi major axis length of the elliptical well
Ry = 50e-10;    % y semi major axis length of the elliptical well
R  = 50e-10;    % radius of the circular well

%% define the mesh
Nx=90;          % Resolution in X axis                 
Ny=70;          % Resoltion in Y axis        
Mx=100e-10;     % Absolute value of the X axis lmit 
My=100e-10;     % Absolute values of the Y axis limit

x=linspace(-Mx,Mx,Nx);  
y=linspace(-My,My,Ny);

[X,Y]=meshgrid(x,y);

%% potential function

% V = rectangular_potential(Wx,Wy,Vo,X,Y);    
% V = elliptical_potential(Wx,Wy,Vo,X,Y);
% V = circular_potential(R,Vo,X,Y);
V = hexagonal_potential(5e-9,Vo,X,Y);

%% Finite Element Method
E1=[];
if length(x)*length(y)>1e4
      N=length(x)*length(y);
      display(strcat('Warning: Take care, H=',num2str(N),'x',num2str(N),'elements'))
end

[E1,psi1] = wave_func_2D(x,y,V,num_sol);

E=nan(num_sol,2);
E(1:length(E1),1)=E1;

%% Display results

disp('Allowed Energies (eV) : ');
display(strcat(num2str(E)));

% visualize the potentials
potential_plotter(x,y,V);

% visualize the wavefunctions
wavefunction_plotter(x,y,V,num_sol,psi1,E);
