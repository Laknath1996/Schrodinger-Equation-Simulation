function V = assymetric_well(z,Wz,Vo)
    V = Vo/2*ones(1,length(z));
    V(-Wz <= z & z <= Wz ) = 0;
    V(z > Wz) = Vo;
end