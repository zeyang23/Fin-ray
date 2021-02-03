%% 单个测试
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
RodA.Newton_conv(TOL);
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

N=20;

theta_series=linspace(theta,0,N);
endpos_series=[linspace(end_pos(1),L0,N);linspace(end_pos(2),0,N)];

Rod= fixedRod(E,L0,wid,thi,n,pos1,[endpos_series(:,N);theta_series(N)]);
Rod.init_exp;

for i=1:N
    theta_last=Rod.conv_theta;
    
    Rod= fixedRod(E,L0,wid,thi,n,pos1,[endpos_series(:,i);theta_series(i)]);
    Rod.init_exp;
    
    Rod.conv_theta=theta_last;
    Rod.update_conv;
    
    TOL=1e-6;
    Rod.Newton_conv(TOL);
    Rod.plot_pos_conv;
end