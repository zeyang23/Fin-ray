% 单个测试
% 21-02-05测试
% 离谱的bug：pos2=[0.3*L0;0.3*L0;60]是可解的，但是对称的下半部分pos2=[0.3*L0;-0.3*L0;-60]却解不出来；
% 求解时L的长度问题，L很容易发散

clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;

pos1=[0*L0;0*L0;0];
pos2=[0.3*L0;0.3*L0;45];
RodA = flexRod(E,L0,wid,thi,n,pos1,pos2);

RodA.init_exp;
TOL=1e-6;
RodA.Newton_conv(TOL);
RodA.plot_pos_conv;


fprintf('\nL')
disp(RodA.Ltotal)