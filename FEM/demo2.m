% 两端固定的柔性板

%%
clear
clc

wid=14e-3;
thi=0.3e-3;
E=200*1e9;

L0=0.255;
n=50;

pdes=[0.155;0.095;pi/3];

R1=planar_nR(E,L0,wid,thi,n,pdes);


f=@(x) cal_balance(x,R1);
options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');
x0=zeros(R1.n_seg+3,1);
[x_solve,fval,exitflag,output] = fsolve(f,x0,options);

R1.theta=x_solve(1:R1.n_seg);
R1.F=x_solve(R1.n_seg+1:R1.n_seg+3);

R1.update;
R1.plot_all;


%%
% 读取ansys输出的形状信息

A=importdata('x1.txt');
x=(A.data(:,2)+A.data(:,5))/1000;

B=importdata('y1.txt');
y=(B.data(:,3)+B.data(:,5))/1000;

hold on

plot(x,y)

%%
function [r,J]=cal_balance(x,Rod)
    Rod.theta=x(1:Rod.n_seg);
    Rod.F=x(Rod.n_seg+1:Rod.n_seg+3);
    
    Rod.update;
    r=Rod.cal_r;
    J=Rod.cal_rdot;
end