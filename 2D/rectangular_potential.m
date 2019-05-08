function V = rectangular_potential(Wx,Wy,Vo,X,Y)
    idx= ( abs(X)< Wx ) .* ( abs(Y)< Wy );
    V = (idx)*0 + (1-idx)*Vo ;
end
