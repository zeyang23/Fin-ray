% 测试 
% 使用fsolve
% 为非线性方程组提供梯度

clear
clc

x0=[0;0];
% options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'Display','off');
options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradients',true);

x = fsolve(@myfun,x0,options);


function [f,grad] = myfun(x)
    f=[4*x(1).^2-20*x(1)+1/4*x(2).^2+8;0.5*x(1).*x(2).^2+2*x(1)-5*x(2)+8];
    if nargout > 1
        grad=[8*x(1)-20 0.5*x(2);0.5*x(2)^2+2 x(1)*x(2)-5];
    end
end