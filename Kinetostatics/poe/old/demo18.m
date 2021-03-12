% 约束下的杆 固定杆长
% 使用牛顿法。计算了梯度
% 使用物体雅可比 Jb

% rotz的输入参数是角度！不是弧度！
% 现在效果就很好了

% 计算连续位姿，实现动画
% 从零状态开始计算连续位姿的方法有问题。目前是倒着算。

clear all;
clc;
addpath(genpath('.'));

%材料的物理参数
I=0.25*pi*1e-12;
E=56*1e9;

%微元数量
n=50;

%杆的长度
L0=1;

delta=L0/n;

%目标位姿
theta=45; % 角度制

N=50;

end_pos=[0.5*L0;0.5*L0];

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

avi_object = VideoWriter('demo18.avi');
avi_object.FrameRate = 10;
open(avi_object);
figure
for I=1:N
    theta_solve=x_all(I,1:n);
    plot_pos2(w_all,q_all,g0,delta,theta_solve,L0)

    M = getframe;
    writeVideo(avi_object,M);

    if I < N
        clf;
    end
end
close(avi_object);

rmpath(genpath('.'));