clear
clc

L0=1;

x0=0;
y0=0.4*L0;

a=1.5*0.12*L0;
b=1.5*0.16*L0;

funx=@(t) a*cos(t)+x0;
funy=@(t) b*sin(t)+y0;

tinterval=[-pi,pi];

hold on
fplot(funx,funy,tinterval)