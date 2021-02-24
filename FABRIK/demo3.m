% 尝试使用FABRIK求解finray

%% 单点情况
clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;

pos1=[0*L0;0*L0;0];
pos2=[0.6*L0;0.6*L0;60];
RodA = fixedRod(E,L0,wid,thi,n,pos1,pos2);

RodA.init_exp;
TOL=1e-6;
RodA.FABRIK_solver(TOL);
RodA.plot_pos_conv;


%% 困难点的单点求解
clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;

pos1=[0*L0;0*L0;0];
pos2=[0.95*L0;0.05*L0;60];
RodA = fixedRod(E,L0,wid,thi,n,pos1,pos2);

RodA.init_exp;
TOL=1e-6;
RodA.FABRIK_solver(TOL);
% RodA.Newton_conv(TOL);
RodA.plot_pos_conv;


%% 直线轨迹测试
clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;

pos1=[0*L0;0*L0;0];

theta=60;
end_pos=[0.6*L0;0.6*L0];

N=5;

theta_series=linspace(0,theta,N);
endpos_series=[linspace(0.7*L0,end_pos(1),N);linspace(0.5*L0,end_pos(2),N)];

Rod= fixedRod(E,L0,wid,thi,n,pos1,[endpos_series(:,1);theta_series(1)]);
Rod.init_exp;
TOL=1e-6;
Rod.Newton_conv(TOL);


for i=1:N
    theta_last=Rod.conv_theta;
    
    Rod= fixedRod(E,L0,wid,thi,n,pos1,[endpos_series(:,i);theta_series(i)]);
    Rod.init_exp;
    
    Rod.conv_theta=theta_last;
    Rod.update_conv;
    
    
    Rod.FABRIK_solver(TOL);
    Rod.plot_pos_conv;
end