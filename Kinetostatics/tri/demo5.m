% 测试，通过3个传感器的数据反解柔性板的形状
% 使用牛顿法求解

% 21-03-22测试结果：效果很离谱，各种奇葩的情况都有
% 可能求解成功；可能求出的形状相近，位置有偏移；可能求解出的形状很离谱；可能求不出来，即不收敛

% 不支持使用zeros(n+3,1)作为初值；
% partial项应该打开

% 21-03-22晚上
% 之前的公式搞错了，应该是sum(theta(p1:q1))，而不是theta(q1)-theta(p1)
% 改回正确的公式后发现效果出奇得好


% 生成3个Delta
clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=40;


% 错误的记录 之前的条件写错了
% pdes=[0.6*L0;0.6*L0;2/3*pi];   % 这组参数的效果非常好
% pdes=[0.6*L0;0.6*L0;1/3*pi];   % 这组参数求出的形状偏得有点远
% pdes=[0.6*L0;0.6*L0;0];        % 这组参数不收敛
% pdes=[0.3*L0;0.3*L0;pi/2];     % 这组参数也不收敛


% 这些参数都是收敛的
% pdes=[0.6*L0;0.6*L0;2/3*pi];   
% pdes=[0.6*L0;0.6*L0;1/3*pi];   
% pdes=[0.6*L0;0.6*L0;0];        
pdes=[0.3*L0;0.3*L0;pi/2];

R1=planar_nR(E,L0,wid,thi,n,pdes);

TOL=1e-6;
R1.Newton(TOL);
R1.plot_all;
hold on

% 传感器的大小和位置
sensor_length=1/8;

p1=fix((1/4-sensor_length/2)*n);
q1=fix((1/4+sensor_length/2)*n);

p2=fix((2/4-sensor_length/2)*n);
q2=fix((2/4+sensor_length/2)*n);

p3=fix((3/4-sensor_length/2)*n);
q3=fix((3/4+sensor_length/2)*n);

Delta1=sum(R1.theta(p1:q1));
Delta2=sum(R1.theta(p2:q2));
Delta3=sum(R1.theta(p3:q3));


% 尝试反解
Rod=planar_nR(E,L0,wid,thi,n,[0;0;0]);

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

% 牛顿法求解N+3方程组
x=zeros(n+3,1);

TOL=1e-6;
k=1;

while(1)
    if k>500
        error("can not converge")
    end
    
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
    
    r=[r1;r2;r3;r4];
    norm(r)
    
    if(norm(r)<TOL)
        fprintf('Newton Method converge: iteration = %d\n',k-1)
        fprintf('norm(e) = %E\n',norm(r))
        break;
    end
    
    J=zeros(n+3);
    
    J(1:n,1:n)=Rod.K_theta-1*Rod.partial;
    J(1:n,n+1:end)=-transpose(Rod.Jacobian);
    J(n+1,1:n)=D1;
    J(n+2,1:n)=D2;
    J(n+3,1:n)=D3;
    
    delta=-pinv(J)*r;
    
    x=x+delta;

    k=k+1;
end

Rod.plot_all