% 调用fsolve求解柔性板形状
% 打开partial项后能够通过Gradient Check

clear
clc

wid=14e-3;
thi=0.3e-3;
E=197*1e9;

L0=0.3;
n=50;

pdes=[0.5*L0;0.5*L0;pi/2];

R1=planar_nR(E,L0,wid,thi,n,pdes);


f=@(x) cal_balance(x,R1);
options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');
x0=zeros(R1.n_seg+3,1);
[x_solve,fval,exitflag,output] = fsolve(f,x0,options);

R1.theta=x_solve(1:R1.n_seg);
R1.F=x_solve(R1.n_seg+1:R1.n_seg+3);

R1.update;
R1.plot_all;


function [r,J]=cal_balance(x,Rod)
    Rod.theta=x(1:Rod.n_seg);
    Rod.F=x(Rod.n_seg+1:Rod.n_seg+3);
    
    Rod.update;
    r=Rod.cal_r;
    J=Rod.cal_rdot;
end