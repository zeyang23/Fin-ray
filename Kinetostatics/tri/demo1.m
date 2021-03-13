% 测试 仅使用三角函数求解两端固定的柔性板
% 21-03-13 放弃指数坐标以后，求解效果非常好

% 21-03-13 总结
% 无论是使用三角函数还是指数坐标，在开启partial(transpose(Jtheta)*F)这项时，都会在自由状态附近产生求解困难。
% 在使用三角函数时，如果关掉partial项，似乎有很强的鲁棒性，能够求解出一些用指数坐标很不好求的位形

clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;

pdes1=[0.2*L0;0.2*L0;pi/3];
pdes2=[0.1*L0;0.1*L0;0]; % 评注：这么离谱的形状都解得出来，三角函数真nb

pdes3=[0.3*L0;0.3*L0;0]; % 评注：这个点解不出来

R1=planar_nR(E,L0,wid,thi,n,pdes1);
% R1.plot_all

% R1.theta=pi/12*ones(n,1);
% 
% R1.cal_pe;
% R1.cal_Jacobian;
% R1.cal_Hessian;

TOL=1e-6;
R1.Newton(TOL);
R1.plot_all;