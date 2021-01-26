function B=Z_adjoint2(xi)
% 01/26更新，原先的Z_adjoint是错的，是用于[omega;v]的
% Z_adjoint2用于[v;omega]
    omega=xi(4:6);
    v=xi(1:3);
    B=zeros(6,6);
    B(1:3,4:6)=xev(omega);
    B(4:6,1:3)=xev(omega);
    B(4:6,4:6)=xev(v);
end