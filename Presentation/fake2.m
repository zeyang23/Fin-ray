% 用来画开题报告里的示意图

clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;


posC1=[-0.8*L0;0*L0;78];
posC2=[-0.5*L0;0.8*L0;20];
RodC = fixedRod(E,0.95*L0,wid,thi,n,posC1,posC2);

RodC.init_exp;
TOL=1e-6;
RodC.Newton_conv(TOL);
RodC.plot_pos;
RodC.cal_pos_all;




posD1=[-0.5*L0;0*L0;102];
posD2=[-0.5*L0;0.8*L0;40];
RodD = fixedRod(E,0.9*L0,wid,thi,n,posD1,posD2);

RodD.init_exp;
TOL=1e-6;
RodD.Newton_conv(TOL);
RodD.plot_pos;
RodD.cal_pos_all;

posA1=[0.2*L0;0*L0;102];
posA2=[-0.1*L0;0.8*L0;160];
RodA = fixedRod(E,0.95*L0,wid,thi,n,posA1,posA2);

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

C0=RodC.pos_all(:,1);
C1=RodC.pos_all(:,8);
C2=RodC.pos_all(:,16);
C3=RodC.pos_all(:,24);
C4=RodC.pos_all(:,32);
C5=RodC.pos_all(:,40);

D0=RodD.pos_all(:,1);
D1=RodD.pos_all(:,8);
D2=RodD.pos_all(:,16);
D3=RodD.pos_all(:,24);
D4=RodD.pos_all(:,32);
D5=RodD.pos_all(:,40);

draw_line(C0,D0)
draw_line(C1,D1)
draw_line(C2,D2)
draw_line(C3,D3)
draw_line(C4,D4)
draw_line(C5,D5)

axis equal
axis([-0.9 0.3 -0.1 1])


object_right=RodB.pos_all(:,20:42);
object_left=RodD.pos_all(:,20:42);

for i =1:22
    A=object_right(:,i);
    A(1)=A(1)-0.01;
    B=object_right(:,i+1);
    B(1)=B(1)-0.01;
    C=object_left(:,i);
    C(1)=C(1)+0.01;
    D=object_left(:,i+1);
    D(1)=D(1)+0.01;
    draw_line2(A,B)
    draw_line2(C,D)
end

R1=object_right(:,end);
R1(1)=R1(1)-0.01;
R2=object_right(:,1);
R2(1)=R2(1)-0.01;
L1=object_left(:,end);
L1(1)=L1(1)+0.01;
L2=object_left(:,1);
L2(1)=L2(1)+0.01;
draw_line2(L1,R1)
draw_line2(L2,R2)

A_end=RodA.pos_all;
B_end=RodB.pos_all;
C_end=RodC.pos_all;
D_end=RodD.pos_all;

save('A_end.mat','A_end')
save('B_end.mat','B_end')
save('C_end.mat','C_end')
save('D_end.mat','D_end')