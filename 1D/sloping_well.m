function V = sloping_well(z,Wz,Vo)
    V = Vo*ones(1,length(z));
    f = @(x) Vo/(4*Wz)*x + Vo/4;
    V(-Wz <= z & z <= Wz ) = f(z(-Wz <= z & z <= Wz));
    V(z > Wz) = Vo;
end