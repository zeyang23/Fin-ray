clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;

pos1=[0.2*L0;0.2*L0;30];

pos2_conv=[0.6*L0;0.6*L0;60];

N=30;
theta1_series=linspace(0,N,30);

for i=1:N
    pos1=[0.2*L0;0.2*L0;theta1_series(i)];
    temp_a=[pos2_conv(1:2);0];
    temp_b=rotz(pos1(3))*temp_a;

    pos2=[temp_b(1:2)+pos1(1:2);pos2_conv(3)+pos1(3)];
    
    RodA = fixedRod(E,L0,wid,thi,n,pos1,pos2);
    RodA.init_exp;
    TOL=1e-6;
    RodA.Newton_conv(TOL);

    RodA.plot_pos;
end