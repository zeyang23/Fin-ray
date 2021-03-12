% 验证横向受力的杆
% 使用牛顿法 计算雅可比矩阵
% F取为定值，可能有问题
% n段微元，尝试取n个旋量，不同于demo8
% F的方向问题，一直很奇怪

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
K=diag((1/c)*ones(n,1));

%外力
F=[0;-20;0;0;0;0];

%初始姿态矩阵
g0=[eye(3),[L;0;0];0,0,0,1];

%所有的w向量构成矩阵
w_all=[];
for i=1:n
    w_all=[w_all,[0;0;1]];
end

%所有的q向量构成矩阵
q_all=zeros(3,n);
for i=1:n
    q_all(:,i)=[1;0;0]*delta*(i-1/2);
end


%求解平衡位置
f=@(theta) cal_balance(K,w_all,q_all,g0,theta,F);
JACOB=@(theta) Jacob_solve(K,w_all,q_all,theta,F);

theta0=linspace(0,pi/100,n);
TOL=1e-6;

theta_solve=Newton_nd(f,JACOB,theta0,TOL);

theta_fsolve=fsolve(f,linspace(0,pi/100,n));

error=norm(theta_solve-theta_fsolve);

clc
disp("difference between Newton and fsolve")
disp(error)

[g_exp,~]=exp_fkine(w_all,q_all,g0,theta_solve);

%得到整个杆的位姿
plot_pos(w_all,q_all,g0,delta,theta_solve)

save('G1.mat','g_exp');
save('Theta1.mat','theta_solve');
save('F1.mat','F');

rmpath(genpath('.'));

function J=Jacob_solve(K,w_all,q_all,theta_all,F)
    L=0.1;
    g0=[eye(3),[L;0;0];0,0,0,1];
    [Jacobs,~,~]=exp_jacob(w_all,q_all,g0,theta_all);
    J=K-partial_J_theta(w_all,q_all,theta_all,Jacobs,F);
end