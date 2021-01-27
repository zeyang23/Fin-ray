% 测试唐学长的数据

% 01/27更新：对求解情况做一个总结
% 雅可比矩阵完全算对也不是件好事
% 在接近原始状态附近时，用正确的雅可比矩阵从theta全0搜是搜不出来的
% KJ直接置0时，对初值的容忍性较高，只不过迭代次数是用正确雅可比时的大概两倍罢了
% 但是正确的雅可比的局限性在于对初值敏感，如果目标姿态接近初始状态，会算不出来
% 需要根据求解的目标姿态来选择是否启用Kj
% 目标姿态如果离初始状态比较远，这时应该启用Kj，好处是迭代次数少，快速收敛
% 目标姿态如果离初始状态近，就应该将Kj置0，如果不置0就算不出来

% 另外，旋量[omega;v]与[v;omega]的问题需要时刻警醒
% 一个Z_adjoint函数写了三版才写对也是服了。。。

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

N=21;

end_pos=[0.6*L0;0.6*L0];

JACOB=@(x) Jacob_constraint_simple2(L0,I,E,n,x);

x_all=interp_solver3(theta,end_pos,N,n,L0,I,E,JACOB);



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