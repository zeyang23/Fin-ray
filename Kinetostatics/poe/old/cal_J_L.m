function J_L=cal_J_L(w_all,q_all,x)
    theta=x(2:end-6);
    n=length(theta)-1;
    
    xi_all=hight(w_all,q_all);
    xi_a=[-cha(w_all,q_all);w_all];
    
    mat_A=zeros(4,4*n+4);
    mat_A(:,1:4)=eye(4);
    H_series=zeros(6,6*n+6);
    H_series(:,1:6)=eye(6);
    
    for i = 2:n+1
        mat_A(:,4*i-3:4*i)= mat_A(:,4*i-7:4*i-4)*mult(xi_all(1:4,4*i-7:4*i-4),theta(i-1));
        H_series(:,6*i-5:6*i)=adjoint(mat_A(:,4*i-3:4*i));
    end
    
    J_L=zeros(6,1);
    for i=1:n+1
        H=H_series(:,6*i-5:6*i);
        Z=Z_adjoint(xi_a(:,i));
        A=theta(i)*eye(6)...
         +0.5*(4-theta(i)*sin(theta(i))-4*cos(theta(i)))*Z...
         +0.5*(4*theta(i)-5*sin(theta(i))+theta(i)*cos(theta(i)))*Z^2 ...
         +0.5*(2-theta(i)*sin(theta(i))-2*cos(theta(i)))*Z^3 ...
         +0.5*(2*theta(i)-3*sin(theta(i))+theta(i)*cos(theta(i)))*Z^4;
        
        if i == n+1
           partial_xi_l=[0;-1;0;0;0;0];
        else
           partial_xi_l=[0;-(i-0.5);0;0;0;0];
        end
        
        J_L=J_L+H*A*partial_xi_l;        
    end
end