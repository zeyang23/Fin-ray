% 21-02-22 更新 新建FABRIK思路
% forward and backward reaching inverse kinematics FABRIK
% 将考虑静力学的模型和纯运动学模型相结合，利用在多机器人课程中学到的瞬时刚性的特点
% 考虑静力学的模型计算复杂度太高，没必要每个点都用考虑静力学的模型来求解逆运动学的问题
% 在某个姿态附近的一个小的邻域内，根据瞬时刚性的概念，即使是冗余机器人，逆运动学也是只有唯一解的

% Extending FABRIK with model constraints 这篇文章对逆运动学的介绍非常棒，将逆运动学的求解方法分为6类

% 灵感：之前对雅可比矩阵求伪逆，暗含着各个关节是平权的，实际上求速度的时候可以解一个带约束的优化问题

clear
clc

n=100;
L=1;

base=[0,0];
target=[0,0.5*L];

d=L/(n-1)*ones(n-1,1);

p0x=transpose(linspace(0,L,n));
p0y=zeros(n,1);

p0=[p0x,p0y];

p=p0;

tol=1e-8;

k=1;

tic;
while(1)
    if norm(p(end,:)-target)<tol
        disp('success')
        disp(k)
        break;
    elseif k>1000
        error('fail to converge')
    else
        p=forward(p,d,target);
        p=backward(p,d,base);
        
        k=k+1;
    end
end
t=toc;
disp(t)

plot_points(p)