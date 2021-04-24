clear
clc

x=-3:1:3;
y=1./(1+x.^2);

u=linspace(-3.2,3.2,100);

v=Spline_natural(x,y,u);

plot(u,v)
hold on
plot(x,y,'*')