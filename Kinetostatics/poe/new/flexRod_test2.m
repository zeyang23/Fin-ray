% 拼接测试
% 中间角度为钝角时难以收敛
clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;

pos1=[0;0;0];

pos2x=0.3*L0;
pos2y=0.3*L0;
pos2theta=60;

% pos2x=0.1*L0;
% pos2y=0.3*L0;
% pos2theta=90;

% pos2x=0.3*L0;
% pos2y=0.1*L0;
% pos2theta=90;


pos2=[pos2x;pos2y;pos2theta];

RodA = flexRod(E,L0,wid,thi,n,pos1,pos2);

RodA.init_exp;
TOL=1e-6;
RodA.Newton_conv(TOL);
RodA.plot_pos;


fprintf('\nL')
disp(RodA.Ltotal)

L_last=L0-RodA.Ltotal;

pos3=[0.6*L0;0.6*L0;90];

RodB=fixedRod(E,L_last,wid,thi,n,pos2,pos3);
RodB.init_exp;
TOL=1e-6;
RodB.Newton_conv(TOL);
RodB.plot_pos;


% 画出相切的两个圆
r=0.01*L0;
vec=r*[cos(pos2theta/180*pi);sin(pos2theta/180*pi);0];
vec1=rotz(90)*vec;
vec2=rotz(-90)*vec;
center1x=vec1(1)+pos2x;
center1y=vec1(2)+pos2y;
center2x=vec2(1)+pos2x;
center2y=vec2(2)+pos2y;

rectangle('Position',[center1x-r,center1y-r,2*r,2*r],'Curvature',[1,1],'linewidth',1,'edgecolor','r')
rectangle('Position',[center2x-r,center2y-r,2*r,2*r],'Curvature',[1,1],'linewidth',1,'edgecolor','r')