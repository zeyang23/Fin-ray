% 测试 
% 使用fmincon
% 含有非线性约束，并且为非线性约束提供梯度

clear
clc

fun = @(x)exp(x(1) + 2*x(2));

x0 = [0 0];
A = []; 
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];

options = optimoptions(@fmincon,'SpecifyConstraintGradient',true);


[x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@ellipseparabola,options);

function [c,ceq,gradc,gradceq] = ellipseparabola(x)
c(1) = x(1)^2/9 + x(2)^2/4 - 1;
c(2) = x(1)^2 - x(2) - 1;
ceq = [];

if nargout > 2
    gradc = [2*x(1)/9, 2*x(1); ...
             x(2)/2, -1];
    gradceq = [];
end
end