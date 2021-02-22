% 拼接测试
% 中间角度为钝角时难以收敛
% 尝试采用初值逼近法
clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;


pos1=[0;0;0];

% 参数
theta_start=45;
theta_end=90;
N=20;
pos2_xy=[0.3*L0;0.3*L0];
pos3=[0.6*L0;0.6*L0;90];



pos2=[pos2_xy;theta_start];
theta_series=linspace(theta_start,theta_end,N);

L_last=L0;
theta_last=zeros(n,1);
F_last=zeros(6,1);


for i = 1:N
    pos2=[pos2_xy;theta_series(i)];
    RodA = flexRod(E,L0,wid,thi,n,pos1,pos2);
    RodA.init_exp;
    
    RodA.Ltotal=L_last;
    RodA.conv_theta=theta_last;
    RodA.conv_F=F_last;
    RodA.update_conv;
    
    TOL=1e-6;
    RodA.Newton_conv(TOL);
    
    L_last=RodA.Ltotal;
    theta_last=RodA.conv_theta;
    F_last=RodA.conv_F;
    
    RodA.plot_pos;
end

figure
RodA = flexRod(E,L0,wid,thi,n,pos1,[pos2_xy;theta_end]);
RodA.init_exp;

RodA.Ltotal=L_last;
RodA.conv_theta=theta_last;
RodA.conv_F=F_last;
RodA.update_conv;

TOL=1e-6;
RodA.Newton_conv(TOL);
RodA.plot_pos;

L_last=L0-RodA.Ltotal;

RodB=fixedRod(E,L_last,wid,thi,n,pos2,pos3);
RodB.init_exp;
TOL=1e-6;
RodB.Newton_conv(TOL);
RodB.plot_pos;