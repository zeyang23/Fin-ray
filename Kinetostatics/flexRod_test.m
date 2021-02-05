% ��������
% 21-02-05����
% ���׵�bug��pos2=[0.3*L0;0.3*L0;60]�ǿɽ�ģ����ǶԳƵ��°벿��pos2=[0.3*L0;-0.3*L0;-60]ȴ�ⲻ������
% ���ʱL�ĳ������⣬L�����׷�ɢ

clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;

pos1=[0*L0;0*L0;0];
pos2=[0.3*L0;-0.3*L0;0];
RodA = flexRod(E,L0,wid,thi,n,pos1,pos2);

RodA.init_exp;
TOL=1e-6;
RodA.Newton_conv(TOL);
RodA.plot_pos_conv;


fprintf('\nL')
disp(RodA.Ltotal)