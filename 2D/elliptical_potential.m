function V = elliptical_potential(Rx,Ry,Vo,X,Y)
    idx= (Y/Ry).^2 + (X/Rx).^2 < 1;
    V = (idx)*0 + (1-idx)*Vo ;
end