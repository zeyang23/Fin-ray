clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;

pdes1=[0.2*L0;0.2*L0;pi/3];
pdes2=[0.1*L0;0.1*L0;0];

R1=planar_nR(E,L0,wid,thi,n,pdes1);

ARRAY=planar_nR.empty;
ARRAY(1)=planar_nR(E,L0,wid,thi,n,pdes1);
ARRAY(2)=planar_nR(E,L0,wid,thi,n,pdes2);

TOL=1e-6;
ARRAY(1).Newton(TOL);
ARRAY(1).plot_all;

hold on

ARRAY(2).Newton(TOL);
ARRAY(2).plot_all;