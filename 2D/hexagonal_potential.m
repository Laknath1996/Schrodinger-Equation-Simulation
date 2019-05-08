function V = hexagonal_potential(A,X,Y,Vo)
    idx1=  (abs(X)<A*sqrt(3)/2);
    idx2=(tan(pi/6)*X+A>Y) .* (tan(pi/6)*X-A<Y) .* (-tan(pi/6)*X-A<Y) .* (-tan(pi/6)*X+A>Y);
    idx=idx1.*idx2;
    V = (idx).*0 + (1-idx).*Vo;
end