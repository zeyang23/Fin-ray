% 约束下的杆，使用内置无梯度优化算法

clear all;
clc;
addpath(genpath('.'));

%材料的物理参数
I=0.25*pi*1e-12;
E=56*1e9;

%微元数量
n=50;

%初始长度
L0=0.1;

%目标位姿
gt=[eye(3),[0.06;0.03;0];0,0,0,1];

f=@(x) cal_constraint(I,E,n,gt,x);

x0=zeros(n+8,1);
x0(1)=L0;
x0(8:end)=linspace(0,0,n+1);

xsolve=fminsearch(f,x0);

theta_final=(xsolve(8:end))';
L_final=xsolve(1);
delta=L_final/n;

%得到整个杆的位姿

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
q_all=[q_all,[1;0;0]*L_final];

pos=[0;0];
g0=[eye(3),[L_final;0;0];0,0,0,1];
[g_exp,~]=exp_fkine(w_all,q_all,g0,theta_final);
for i =1:n
    g0i=[eye(3),[(i-0.5)*delta;0;0];0,0,0,1];
    [gsti,~]=exp_fkine(w_all(:,1:i),q_all(:,1:i),g0i,theta_final(:,1:i));
    pos=[pos,gsti(1:2,4)];
end
pos=[pos,g_exp(1:2,4)];

plot(pos(1,:),pos(2,:))

rmpath(genpath('.'));