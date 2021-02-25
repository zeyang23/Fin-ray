% 尝试使用FABRIK生成合适的初值，以提高牛顿法的收敛速度，以及减少在给定初值下牛顿法不收敛的问题

clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;

pos1=[0*L0;0*L0;0];
pos2=[0.6*L0;-0.6*L0;-60];
RodA = fixedRod(E,L0,wid,thi,n,pos1,pos2);

RodA.init_exp;
TOL=1e-6;
RodA.FABRIK_solver(TOL);
RodA.update_conv;
RodA.Newton_conv(TOL);

RodA.plot_pos_conv;