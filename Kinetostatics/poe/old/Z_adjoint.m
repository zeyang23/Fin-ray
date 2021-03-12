function B=Z_adjoint(xi)
    omega=xi(4:6);
    v=xi(1:3);
    B=zeros(6,6);
    B(1:3,1:3)=xev(omega);
    B(4:6,4:6)=xev(omega);
    B(4:6,1:3)=xev(v);
end