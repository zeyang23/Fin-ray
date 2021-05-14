%
% 长度单位取mm
% 杨氏模量单位取Mpa

% 注意，算出的力的单位为kN，力矩的单位为Nm


clear
clc


wid=14;
thi=0.2;
E=197;

L0=350;
n=100;


pdes=[L0-0.4*L0;0.4*L0;pi/3];
x0=zeros(n+3,1);


Rod=planar_nR(E,L0,wid,thi,n,pdes);


f=@(x) cal_balance(x,Rod);
options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

[x_solve,fval,exitflag,output] = fsolve(f,x0,options);

Rod.theta=x_solve(1:Rod.n_seg);
Rod.F=x_solve(Rod.n_seg+1:Rod.n_seg+3);

Rod.update;
Rod.plot_all;

disp(Rod.F);


function [r,J]=cal_balance(x,Rod)
    Rod.theta=x(1:Rod.n_seg);
    Rod.F=x(Rod.n_seg+1:Rod.n_seg+3);
    
    Rod.update;
    r=Rod.cal_r;
    J=Rod.cal_rdot;
end