function C=solve_triangle2(A,B,psi)
    AB=B-A;
    L=norm(AB);
    theta=atan2(AB(2),AB(1));
    phi=psi-theta;
    C=[B(1)-L*cos(phi),B(2)+L*sin(phi)];
end
    