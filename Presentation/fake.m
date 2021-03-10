% 用来画开题报告里的示意图

clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;


pos1=[0.2*L0;0*L0;102];
pos2=[-0.1*L0;0.8*L0;160];
RodA = fixedRod(E,0.95*L0,wid,thi,n,pos1,pos2);

RodA.init_exp;
TOL=1e-6;
RodA.Newton_conv(TOL);
RodA.plot_pos;
RodA.cal_pos_all;




posB1=[-0.1*L0;0*L0;78];
posB2=[-0.1*L0;0.8*L0;140];
RodB = fixedRod(E,0.9*L0,wid,thi,n,posB1,posB2);

RodB.init_exp;
TOL=1e-6;
RodB.Newton_conv(TOL);
RodB.plot_pos;
RodB.cal_pos_all;

axis equal
axis([-0.5 0.5 -0.1 1])

xticks(linspace(-0.5,0.5,11))
yticks(linspace(0,1,11))


x=-0.07;
y=0.5;
r=0.15;
rectangle('Position',[x-r,y-r,2*r,2*r],'Curvature',[1,1],'linewidth',1,'edgecolor','r')

A0=RodA.pos_all(:,1);
A1=RodA.pos_all(:,8);
A2=RodA.pos_all(:,16);
A3=RodA.pos_all(:,24);
A4=RodA.pos_all(:,32);
A5=RodA.pos_all(:,40);

B0=RodB.pos_all(:,1);
B1=RodB.pos_all(:,8);
B2=RodB.pos_all(:,16);
B3=RodB.pos_all(:,24);
B4=RodB.pos_all(:,32);
B5=RodB.pos_all(:,40);

draw_line(A0,B0)
draw_line(A1,B1)
draw_line(A2,B2)
draw_line(A3,B3)
draw_line(A4,B4)
draw_line(A5,B5)