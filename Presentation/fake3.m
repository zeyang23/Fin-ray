% 用来画开题报告里的示意图

clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;

trans=0.11*L0;


THETA=81.3;
THETA2=180-THETA;

posA1=[trans+0.2*L0;0*L0;THETA2];

alpha=THETA2/180*pi;

posA2=[trans+0.2*L0+L0*cos(alpha);L0*sin(alpha);THETA2];
RodA = fixedRod(E,L0,wid,thi,n,posA1,posA2);

RodA.init_exp;
TOL=1e-6;
RodA.Newton_conv(TOL);
RodA.plot_pos;
RodA.cal_pos_all;

posB1=[trans-0.1*L0;0*L0;THETA];
beta=THETA/180*pi;
posB2=[trans-0.1*L0+L0*cos(beta);L0*sin(beta);THETA];
RodB = fixedRod(E,L0,wid,thi,n,posB1,posB2);

RodB.init_exp;
TOL=1e-6;
RodB.Newton_conv(TOL);
RodB.plot_pos;
RodB.cal_pos_all;


trans2=-0.92*L0;

posC1=[trans2+trans+0.2*L0;0*L0;THETA2];

alpha=THETA2/180*pi;

posC2=[trans2+trans+0.2*L0+L0*cos(alpha);L0*sin(alpha);THETA2];
RodC = fixedRod(E,L0,wid,thi,n,posC1,posC2);


RodC.init_exp;
TOL=1e-6;
RodC.Newton_conv(TOL);
RodC.plot_pos;
RodC.cal_pos_all;

posD1=[trans2+trans-0.1*L0;0*L0;THETA];
beta=THETA/180*pi;
posD2=[trans2+trans-0.1*L0+L0*cos(beta);L0*sin(beta);THETA];
RodD = fixedRod(E,L0,wid,thi,n,posD1,posD2);

RodD.init_exp;
TOL=1e-6;
RodD.Newton_conv(TOL);
RodD.plot_pos;
RodD.cal_pos_all;

A_start=RodA.pos_all;
B_start=RodB.pos_all;
C_start=RodC.pos_all;
D_start=RodD.pos_all;

save('A_start.mat','A_start')
save('B_start.mat','B_start')
save('C_start.mat','C_start')
save('D_start.mat','D_start')