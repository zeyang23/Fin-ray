% 约束下的杆 固定杆长
% 使用牛顿法。计算了梯度
% 使用物体雅可比 Jb

% rotz的输入参数是角度！不是弧度！
% 现在效果就很好了

% 尝试使用空间雅可比

% 发现问题：最后一个关节的求解是错的

clear all;
clc;
addpath(genpath('.'));

%材料的物理参数
wid=5e-3;
thi=1e-3;
I=1/12*thi.^3*wid;
E=197*1e9;

%微元数量
n=20;

%杆的长度
L0=0.2;

delta=L0/n;

%目标位姿
theta=60; % 角度制

%目标位姿
gt=[rotz(theta),[0.6*L0;0.6*L0;0];0,0,0,1];

f=@(x) cal_constraint_simple(L0,I,E,n,gt,x);

x0=zeros(n+6,1);

x0(1:end-6)=linspace(0,0,n);

x0(end-5:end)=0*ones(6,1);

% load('x_demo13.mat')
% x0=xsolve;

JACOB=@(x) Jacob_constraint_simple(L0,I,E,n,x);

TOL=1e-3;

% x_fsolve=fsolve(f,x0);

xsolve=Newton_nd(f,JACOB,x0,TOL);

theta_solve=(xsolve(1:end-6))';


%初始姿态矩阵
g0=[eye(3),[L0;0;0];0,0,0,1];

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

[g_exp,~]=exp_fkine(w_all,q_all,g0,theta_solve);

%得到整个杆的位姿
plot_pos(w_all,q_all,g0,delta,theta_solve)

rmpath(genpath('.'));