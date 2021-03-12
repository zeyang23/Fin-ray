clear
clc

f=@(x) [4*x(1).^2-20*x(1)+1/4*x(2).^2+8;0.5*x(1).*x(2).^2+2*x(1)-5*x(2)+8];
Jacob=@(x) [8*x(1)-20 0.5*x(2);
            0.5*x(2)^2+2 x(1)*x(2)-5];

x0=[0;0];
TOL=1e-6;

x=Newton_nd(f,Jacob,x0,TOL)

fsolve(f,x0)