function A=partial_J_theta(w_all,q_all,theta_all,Jacobs,F)
    xi_a=[-cha(w_all,q_all);w_all];
    n=size(xi_a,2);
    A=zeros(n,n);
    for J=2:n
        for K=1:J-1
            A(J,K)=transpose(Z_adjoint(Jacobs(:,K))*Jacobs(:,J))*F;
        end
    end
end