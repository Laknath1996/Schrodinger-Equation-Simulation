function[E,psi]=wave_func_2D(x,y,V,n)

% constants
h_bar = 1.055e-34;
m_e = 9.11e-31;
e = 1.602e-19;

Nx=length(x);
Ny=length(y);

dx=x(2)-x(1);
dy=y(2)-y(1);

Axy = ones(1,(Nx-1)*Ny);
DX2 = (-2)*diag(ones(1,Ny*Nx)) + (1)*diag(Axy,-Ny) + (1)*diag(Axy,Ny);

AA=ones(1,Nx*Ny);
BB=ones(1,Nx*Ny-1);
BB(Ny:Ny:end)=0;

DY2=(-2)*diag(AA) + (1)*diag(BB,-1) + (1)*diag(BB,1);

H=(-h_bar^2/(2*m_e)) * ( DX2/dx^2 + DY2/dy^2 ) +  diag(V(:)*e) ;

H=sparse(H);
[PSI,Energy] = eigs(H,n,'SM');
E = diag(Energy)/e;

for i=1:n
    psi_temp=reshape(PSI(:,i),Ny,Nx);
    psi(:,:,i) = psi_temp / sqrt( trapz( y' , trapz(x,abs(psi_temp).^2 ,2) , 1 )  );  % normalisation of the wave function psi
end

psi=psi(:,:,end:-1:1);
E=E(end:-1:1);

end
