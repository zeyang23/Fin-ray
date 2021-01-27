function B=Z_adjoint2(xi)
% 01/26���£�ԭ�ȵ�Z_adjoint�Ǵ�ģ�������[omega;v]��
% Z_adjoint2����[v;omega]
% 01/26����,Z_adjoint2Ҳ�Ǵ�ģ�û�к�F��Ӧ��ԭ���е�FΪ[m;f]�����õ���[f;m]
    omega=xi(4:6);
    v=xi(1:3);
    B=zeros(6,6);
    B(1:3,1:3)=xev(omega);
    B(4:6,4:6)=xev(omega);
    B(1:3,4:6)=xev(v);
end