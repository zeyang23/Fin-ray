% 测试唐学长的数据

clear all;
clc;
addpath(genpath('.'));

%材料的物理参数
wid=5e-3;
thi=1e-3;
I=1/12*thi.^3*wid;
E=197*1e9;

%微元数量
n=50;

%杆的长度
L0=0.2;

delta=L0/n;

%目标位姿
theta=60; % 角度制

N=20;

end_pos=[0.6*L0;0.6*L0];

JACOB=@(x) Jacob_constraint_simple2(L0,I,E,n,x);

x_all=interp_solver3(theta,end_pos,N,n,L0,I,E,JACOB);

JACOB=@(x) Jacob_constraint_simple2(L0,I,E,n,x);


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


figure
for i=1:N
    theta_solve=x_all(i,1:n);
    plot_pos2(w_all,q_all,g0,delta,theta_solve,L0)

end

rmpath(genpath('.'));