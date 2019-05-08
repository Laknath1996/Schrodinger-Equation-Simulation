function V = double_well(z,Wz,Vo)
    V = Vo/2*ones(1,length(z));
    V(-Wz <= z & z <= Wz ) = Vo;
    V(-2*Wz <= z & z < -Wz) = 0;
    V(Wz <= z & z < 2*Wz) = 0;
end
