% 验证横向受力的杆
% 翻转雅可比矩阵，demo4中的J未经翻转，可能是错的
% demo7与demo4得到的结果相同
% y方向的正负没搞明白

clear all;
clc;
addpath(genpath('.'));

%离散的微元个数
n=50;

%柔性杆的总长度
L=0.1;

%微元的长度
delta=L/n;

%材料的物理参数
I=0.25*pi*1e-12;
E=56*1e9;
c=delta/(E*I);

%刚度矩阵
K=diag((1/c)*ones(n+1,1));

%外力
F=[0;0;0;0;-10;0];

%初始姿态矩阵
g0=[eye(3),[L;0;0];0,0,0,1];

%所有的w向量构成矩阵
w_all=[];
for i=1:n
    w_all=[w_all,[0;0;1]];
end
w_all=[w_all,[0;0;1]];

%所有的q向量构成矩阵
q_all=zeros(3,n);
for i=1:n
    q_all(:,i)=[1;0;0]*delta*(i-1/2);
end
q_all=[q_all,[1;0;0]*L];

%求解平衡位置
f=@(theta) cal_balance2(K,w_all,q_all,g0,theta,F);
theta_solve=fsolve(f,linspace(0,pi/100,n+1));

[g_exp,~]=exp_fkine(w_all,q_all,g0,theta_solve);

%得到整个杆的位姿
pos=[0;0];
for i =1:n
    g0i=[eye(3),[(i-0.5)*delta;0;0];0,0,0,1];
    [gsti,~]=exp_fkine(w_all(:,1:i),q_all(:,1:i),g0i,theta_solve(:,1:i));
    pos=[pos,gsti(1:2,4)];
end
pos=[pos,g_exp(1:2,4)];

plot(pos(1,:),pos(2,:),'-o','MarkerSize',2)

rmpath(genpath('.'));