function V = potential_func(z,Wz,f,sel)
    if sel == 0
        V = f(z);
    end
    if sel == 1
        V = zeros(1,length(z));
        V(z <= -Wz) = f(z(z < -Wz));
        V(z >= Wz) = f(-z(z > Wz));
    end
    if sel == 2
        V = f(Wz)*ones(1,length(z));
        V(z <= Wz & z >= -Wz) = f(z(z <= Wz & z >= -Wz));
    end
end
