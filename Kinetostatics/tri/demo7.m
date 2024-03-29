% 测试，通过多于3个传感器的数据反解柔性板的形状
% 使用fsolve求解最小二乘问题，使用的其实是Levenberg-Marquardt algorithm，和nonlsq中是一样的

% 21-03-22 4个传感器在某些情况下是能求解的
% 不提供梯度是可以解的，提供梯度后收敛速度显著变快。注意，这里提供的梯度要打开partial项，否则CheckGradient会报错
% 存在初值敏感性。存在解出来的值相对原始值大幅偏离的可能性。

% 21-03-22晚上
% 之前的公式搞错了，应该是sum(theta(p1:q1))，而不是theta(q1)-theta(p1)
% 改回正确的公式后发现效果出奇得好


% 生成4个Delta
clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=40;


% 这些参数都是收敛的
pdes=[0.6*L0;0.6*L0;2/3*pi];   
% pdes=[0.6*L0;0.6*L0;1/3*pi];   
% pdes=[0.6*L0;0.6*L0;0];        
% pdes=[0.3*L0;0.3*L0;pi/2];

R1=planar_nR(E,L0,wid,thi,n,pdes);

TOL=1e-6;
R1.Newton(TOL);
R1.plot_all;
hold on

% 传感器的大小和位置
sensor_length=1/8;

p1=fix((1/5-sensor_length/2)*n);
q1=fix((1/5+sensor_length/2)*n);

p2=fix((2/5-sensor_length/2)*n);
q2=fix((2/5+sensor_length/2)*n);

p3=fix((3/5-sensor_length/2)*n);
q3=fix((3/5+sensor_length/2)*n);

p4=fix((4/5-sensor_length/2)*n);
q4=fix((4/5+sensor_length/2)*n);

Delta1=sum(R1.theta(p1:q1));
Delta2=sum(R1.theta(p2:q2));
Delta3=sum(R1.theta(p3:q3));
Delta4=sum(R1.theta(p4:q4));

D1=zeros(1,n);
for i=p1:q1
    D1(i)=1;
end

D2=zeros(1,n);
for i=p2:q2
    D2(i)=1;
end

D3=zeros(1,n);
for i=p3:q3
    D3(i)=1;
end

D4=zeros(1,n);
for i=p4:q4
    D4(i)=1;
end


% 使用fsolve求解

x0=zeros(n+3,1);
options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',true);

f=@(x) myfun(x,E,L0,wid,thi,n,D1,D2,D3,D4,Delta1,Delta2,Delta3,Delta4);
[x,fval,exitflag,output] = fsolve(f,x0,options);
% [x,fval,exitflag,output]= fsolve(f,x0);

theta_real=R1.theta;
x_real=[R1.theta;R1.F];

disp(norm(x(1:n)-theta_real));
disp(norm(x-x_real));

Rod_solve=planar_nR(E,L0,wid,thi,n,pdes);
Rod_solve.theta=x(1:n);
Rod_solve.F=x(end-2:end);

Rod_solve.update;
Rod_solve.plot_all;

function [f,grad] = myfun(x,E,L0,wid,thi,n,D1,D2,D3,D4,Delta1,Delta2,Delta3,Delta4)
    Rod=planar_nR(E,L0,wid,thi,n,[0;0;0]);
    
    theta=zeros(n,1);
    F=zeros(3,1);
    
    theta=x(1:n);
    F=x(end-2:end);
    
    Rod.theta=theta;
    Rod.F=F;
    
    Rod.update;
    
    r1=Rod.K_theta*theta-transpose(Rod.Jacobian)*F;
    r2=D1*theta-Delta1;
    r3=D2*theta-Delta2;
    r4=D3*theta-Delta3;
    r5=D4*theta-Delta4;
    
    f=[r1;r2;r3;r4;r5];
    
    if nargout > 1
        J=zeros(n+3);
    
        J(1:n,1:n)=Rod.K_theta-1*Rod.partial;
        J(1:n,n+1:end)=-transpose(Rod.Jacobian);
        J(n+1,1:n)=D1;
        J(n+2,1:n)=D2;
        J(n+3,1:n)=D3;
        J(n+4,1:n)=D4;
        grad=J;
    end
end