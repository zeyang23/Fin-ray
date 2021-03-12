function [xi_all] = hight(w_all,q_all)       %%将R^6→se(3)，并补上垂直叉乘部分
    xi_all=zeros(4,4*size(w_all,2));
    for a=1:size(w_all,2)
        xi_all(1:3,4*(a-1)+1:4*(a-1)+3)=[0,-w_all(3,a),w_all(2,a);w_all(3,a),0,-w_all(1,a);-w_all(2,a),w_all(1,a),0];
        xi_all(1:3,4*a)=-cha(w_all(1:3,a),q_all(1:3,a));
    end
end