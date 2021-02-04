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

theta=45;
pos2=[0.2*L0;0.2*L0;theta];

N=20;
theta_series=linspace(theta,120,N);

L_last=L0;
theta_last=zeros(n,1);
F_last=zeros(6,1);


for i = 1:N
    pos2=[0.2*L0;0.2*L0;theta_series(i)];
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
    
%     RodA.plot_pos;
end
%%
RodA = flexRod(E,L0,wid,thi,n,pos1,[0.2*L0;0.2*L0;120]);
RodA.init_exp;

RodA.Ltotal=L_last;
RodA.conv_theta=theta_last;
RodA.conv_F=F_last;
RodA.update_conv;

TOL=1e-6;
RodA.Newton_conv(TOL);
RodA.plot_pos;

L_last=L0-RodA.Ltotal;

pos3=[0.6;0.6;90];

RodB=fixedRod(E,L_last,wid,thi,n,pos2,pos3);
RodB.init_exp;
TOL=1e-6;
RodB.Newton_conv(TOL);
RodB.plot_pos;
