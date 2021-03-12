function [gi_mult] = mult(xi_all,theta_all)
    gi_mult=eye(4);
    for a=1:length(theta_all)
        w=xi_all(1:3,4*(a-1)+1:4*(a-1)+3);
        v=xi_all(1:3,4*a);
        exp_w=(eye(3)+w*sin(theta_all(a))+w*w*(1-cos(theta_all(a))));
        x=(eye(3)-exp_w)*w*v+[w(3,2),w(1,3),w(2,1)]*v*transpose([w(3,2),w(1,3),w(2,1)])*theta_all(a);
        gi_mult=gi_mult*[exp_w,x;0 0 0 1];
    end
end