function V = circular_potential(R,Vo,X,Y)
    idx= (Y/R).^2 + (X/R).^2 < 1;
    V = (idx)*0 + (1-idx)*Vo ;
end