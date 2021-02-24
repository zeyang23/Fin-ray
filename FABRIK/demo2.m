% ²âÊÔ·â×°µÄº¯Êý simpleFABRIK

clear
clc

n=100;
L=1;

base=[0,0];
target=[0.6*L,0.5*L];

d=L/(n-1)*ones(n-1,1);

p0x=transpose(linspace(0,L,n));
p0y=zeros(n,1);

p0=[p0x,p0y];

tol=1e-8;

p_solve=simpFABRIK(base,target,d,p0,tol);

plot_points(p_solve)