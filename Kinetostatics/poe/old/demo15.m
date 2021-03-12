% 约束下的杆 固定杆长
% 使用牛顿法。计算了梯度
% 使用受力下的杆的结果来校验雅可比的正确性 来自demo10
% 如果使用物体雅可比，效果比空间雅可比好很多
% 问题：如果自己随便指定末端姿态和位置，还是算不出来，必须得是从受力的杆生成的真实结果
% 尝试不指定末端位姿，仅指定末端位置

clear all;
clc;
addpath(genpath('.'));

%材料的物理参数
I=0.25*pi*1e-12;
E=56*1e9;

%微元数量
n=50;

%杆的长度
L0=0.1;

delta=L0/n;

gt=[rotz(45/180*pi),[0.6*L0;0.8*L0;0];0,0,0,1];

f=@(x) cal_constraint_simple3(L0,I,E,n,gt,x);

x0=zeros(n+6,1);


JACOB=@(x) Jacob_constraint_simple2(L0,I,E,n,x);

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

[g_solve,~]=exp_fkine(w_all,q_all,g0,theta_solve);
disp('g_solve')
disp(g_solve)

%得到整个杆的位姿
plot_pos(w_all,q_all,g0,delta,theta_solve)

rmpath(genpath('.'));