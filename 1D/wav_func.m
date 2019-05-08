function psi = wav_func(E,V,dz,N)
    h_bar = 1.055e-34;
    m_e = 9.11e-31;
    e = 1.602e-19;

    psi = zeros(1,N);
    psi(2) = 1;
    
    for n = 2:N-1
       kc = (2*m_e*e/h_bar^2)*(E - V(n))*dz^2;
       psi(n+1) = (2 - kc)*psi(n) - psi(n-1);
    end
    
    Area = trapz(psi.*psi);
    psi = psi /sqrt(Area);
    
end
