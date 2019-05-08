function V = stepped_well(z,Wz,Vo)
    V = Vo*ones(1,length(z));
    V(-Wz <= z & z <= Wz ) = 0;
    V(z >= 0 & z <= Wz ) = Vo/2;
end

    