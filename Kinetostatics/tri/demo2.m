% 直线轨迹测试

% 测试 仅使用三角函数求解两端固定的柔性板
% 21-03-13 放弃指数坐标以后，求解效果非常好

clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;

% pdes=[0.2*L0;0.2*L0;pi/3];

theta=pi/3;

pdes=[0.6*L0;0.6*L0;theta];

R1=planar_nR(E,L0,wid,thi,n,pdes);

N=20;

theta_series=linspace(theta,0,N);
endpos_series=[linspace(pdes(1),L0,N);linspace(pdes(2),0,N)];

Rod= planar_nR(E,L0,wid,thi,n,[endpos_series(:,N);theta_series(N)]);

TOL=1e-6;
R1.Newton(TOL);

for i=1:N
    theta_last=Rod.theta;
    
    Rod= planar_nR(E,L0,wid,thi,n,[endpos_series(:,i);theta_series(i)]);
    
%     Rod.theta=theta_last;
%     Rod.update;
    
    TOL=1e-6;
    Rod.Newton(TOL);
    Rod.plot_all;
    hold on
end