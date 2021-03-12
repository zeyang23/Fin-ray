% 对比测试

%% 下半部分作为flex 上半部分作为fix
clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;

pos1=[0;0;0];
pos2=[0.3;0.3;60];
RodA = flexRod(E,L0,wid,thi,n,pos1,pos2);

RodA.init_exp;
TOL=1e-6;
RodA.Newton_conv(TOL);
RodA.plot_pos;


fprintf('\nL')
disp(RodA.Ltotal)

L_last=L0-RodA.Ltotal;

pos3=[0.6;0.6;90];

RodB=fixedRod(E,L_last,wid,thi,n,pos2,pos3);
RodB.init_exp;
TOL=1e-6;
RodB.Newton_conv(TOL);
RodB.plot_pos;

%% 上半部分作为flex 下半部分作为fix
clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;

pos1=[0.6;0.6;90];
pos3=[0;0;0];
pos2=[0.3;0.3;-120];

RodA = flexRod(E,L0,wid,thi,n,pos1,pos2);

RodA.init_exp;

RodA.update_conv;

TOL=1e-6;
RodA.Newton_conv(TOL);
RodA.plot_pos;


fprintf('\nL')
disp(RodA.Ltotal)

L_last=L0-RodA.Ltotal;



RodB=fixedRod(E,L_last,wid,thi,n,pos2,pos3);
RodB.init_exp;
TOL=1e-6;
RodB.Newton_conv(TOL);
RodB.plot_pos;