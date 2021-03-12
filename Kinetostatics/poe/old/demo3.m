%验证杆的情况

clear all;
clc;
addpath(genpath('.'));

%离散的微元个数
n=50;

%柔性杆的总长度
L=1;

%微元的长度
delta=L/n;

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

theta_all=pi/100*ones(1,n);
theta_all=[theta_all,0];

[g_exp,~]=exp_fkine(w_all,q_all,g0,theta_all);


%得到整个杆的位姿
pos=[0;0];
for i =1:n
    g0i=[eye(3),[(i-0.5)*delta;0;0];0,0,0,1];
    [gsti,~]=exp_fkine(w_all(:,1:i),q_all(:,1:i),g0i,theta_all(:,1:i));
    pos=[pos,gsti(1:2,4)];
end

pos=[pos,g_exp(1:2,4)];

plot(pos(1,:),pos(2,:),'-o','MarkerSize',2)

rmpath(genpath('.'));