function [jacob0,jacobe,gap] = exp_jacob(w_all,q_all,g0,theta_all)
    n=length(theta_all);
    jacob0=zeros(6,size(w_all,2));
    xi_all=hight(w_all,q_all);                          %½«R^6¡úse(3)
    xi_a=[-cha(w_all,q_all);w_all];
    
    mat_A=zeros(4,4*n);
    mat_A(:,1:4)=eye(4);
    for i = 2:n
        mat_A(:,4*i-3:4*i)= mat_A(:,4*i-7:4*i-4)*mult(xi_all(1:4,4*i-7:4*i-4),theta_all(i-1));
    end
    
    A_end=mat_A(:,4*n-3:4*n)*mult(xi_all(:,4*n-3:4*n),theta_all(n))*g0;
    
    for a=1:size(w_all,2)
        jacob0(:,a)=adjoint(mat_A(:,4*(a-1)+1:4*a))*xi_a(:,a);
    end
    
    jacobe=adjoint(A_end)\jacob0;
    gap=gap_mat(A_end);
end