% author : Ashwin de Silva - 150103P

%% Solutions to Three Dimensional Schrödinger Equation - H Atom

close all;
clear all;

%% user inputs

r_max = input('max radial distance (default 10e-10 m), r_max =  ');     % maximum radial distance considered
disp('   ');
L = input('orbital quantum number (default 0), L =  ');                 % orbital quantum number
disp('   ');
m_L = input('magnetic quantum number (default 0), m_L =  ');            % magnetic quantum number
disp('   ');


num = 1201; % resolution in the r axis
Z = 1;      % nuclear charge

%% constants

hbar = 1.055e-34;      % J.s
e = 1.602e-19;         % C
me = 9.109e-31;        % kg
eps0 = 8.854e-12;      % F/m

%% potential well

% potential energy in electron volts (eV)
r_min = 1e-15;    % default 1e-15
r = linspace(r_min,r_max, num);
dr = r(2)-r(1);
dr2 = dr^2;

% Coulomb term
K = -Z*e/(4*pi*eps0);
U_c = K./r;

% Angular momnetum term
U_L = (hbar^2*L*(L+1)/(2*me*e))./r.^2;

% Effective potential energy
U = U_c + U_L;

U(U < -2000) = -2000;     % minimum value that U can take 
U(U > 1000) = 1000;       % maximum value that U can take

for cn =1:(num-2)
    U_matrix(cn,cn) = U(cn+1);
end

%% Solve Schrodinger Equation for Radial Axis

off = ones(num-3,1);                 
SD_matrix = (-2*eye(num-2) + diag(off,1) + diag(off,-1))/dr2;   % second derivative matrix

% Make KE Matrix
K_matrix = -hbar^2/(2*me*e) * SD_matrix;                        % Kinetic Energy Matrix

% Make Hamiltonian Matrix
H_matrix = K_matrix + U_matrix;                                 % Hamlitonian (H)

% Find Eignevalues E_n and Eigenfunctions psi_N 
[e_funct e_values] = eig(H_matrix);                             % Find the eigenvalues and eigenvectors of the H

% All Eigenvalues 1, 2 , ... n  where E_N < 0
flag = 0;
n = 1;
while flag == 0
    E(n) = e_values(n,n);
    if E(n) > 0
        flag = 1; 
    end
    n = n + 1;
end
E(n-1) = [];
n = n-2;

% Get the eigenvectors
for cn = 1 : n
psi(:,cn) = [0; e_funct(:,cn); 0];
area = simpson1d((psi(:,cn).* psi(:,cn))',r_min,r_max);
psi(:,cn) = psi(:,cn)/sqrt(area);
prob(:,cn) = psi(:,cn) .* psi(:,cn);
end

disp('  ')
disp('Enter Principal Quantum Number for calculation of')
disp('   expectation values and graphical display')
pqn = input('Enter Principal Quantum Number (n > L),  n  =  ');
qn = pqn-L;

Psi = psi(:,qn)';               % radial wavefunction as specified by the pqn
Prob_density = Psi .* Psi;      % radial probability density

%% generate Probability Clouds

nump = 201;         % resolution
maxp = r_max;       

rpd = Prob_density;
EB = -E(pqn-L);     % associated energy of the state

dp = 2*maxp/(nump-1);
xp = -maxp : dp : maxp; 
yp = -maxp : dp : maxp;

k = 0;

for c1 = 1 : nump
    for c2 = 1 : nump
   
        rp = sqrt(xp(c1)^2 + yp(c2)^2);
   
        if rp < maxp   
            % find the theta value for the given xp and yp
            if xp(c1) >= 0 & yp(c2) >= 0
               theta = atan(abs(yp(c2)/(xp(c1)+eps)));
            end
            if xp(c1) <= 0 & yp(c2) >= 0
               theta = pi - atan(abs(yp(c2)/(xp(c1)+eps)));
            end
            if xp(c1) <= 0 & yp(c2) <= 0
               theta = pi + atan(abs(yp(c2)/(xp(c1)+eps)));
            end
            if xp(c1) >= 0 & yp(c2) <= 0
               theta = 2*pi - atan(abs(yp(c2)/(xp(c1)+eps)));
            end

            leg01 = legendre(L,cos(theta));     % get the Lth order legendre polynomials for the theta value (L+1 polynomials)
            leg02 = leg01(m_L+1,:);             % choose the polynomial from the value vector according to the mL

            N = round(1 + (num-1)*(rp/r_max));  % get the mapped index to the radial wavefunction
            probcloud(c1,c2) = rpd(N)*leg02^2;  % find the prob in finding an electron in (xp,yp)
        else
            probcloud(c1,c2) = 0;

        end 

    end
end

probcloud = probcloud./max(max(probcloud));     % normalize the porbablity cloud
probS = probcloud .^0.5;

%% Generate Plots

s = sprintf('Electron density:  n = %.0g    l = %.0g   m = %.0g     EB = %0.4g  eV',pqn,L,m_L,EB);

figure(99)
set(gcf,'Units','normalized','Position',[0.1 0.2 0.8 0.5]);
fs = 16;
set(gca,'fontsize',fs);
subplot(1,3,1)
xx = xp .* 1e10;  yy = xp .*1e10;
pcolor(xx,yy,probS);
colormap jet
colorbar off;
%shading flat;
shading interp
axis square
%title(s);
xlabel('x  (angstroms)  ','fontsize',fs);
ylabel('y   (angstroms) ','fontsize',fs);
%axis off
set(gca,'fontsize',fs);

subplot(1,3,2)
surf(xp,yp,probcloud);
colorbar off;
shading interp;
s = sprintf('Electron density:  n = %.0g    L = %.0g   mL = %.0g     EB = %0.4g  eV',pqn,L,m_L,EB);
title(s,'fontsize',fs);
%xlabel('x  (ao)  ');
%ylabel('y   (ao) ');
axis off
view(-83,18)
hold off;

subplot(1,3,3);
point_cloud(probcloud);



