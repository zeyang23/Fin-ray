% 倾斜 两端固定的柔性板
% 调用fsolve求解

clear
clc

wid=14e-3;
thi=0.3e-3;
E=197*1e9;

L0=0.3;
n=50;

xA=0.1*L0;
yA=0.1*L0;
alpha=pi/3;

pA=[xA;yA;alpha];
pdes=[0.5*L0;0.5*L0;pi/2];

R1=planar_nR(E,L0,wid,thi,n,pdes);

f=@(x) cal_balance(x,R1,pA,pdes);
options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');
x0=zeros(R1.n_seg+3,1);
[x_solve,fval,exitflag,output] = fsolve(f,x0,options);

R1.theta=x_solve(1:R1.n_seg);
R1.F=x_solve(R1.n_seg+1:R1.n_seg+3);

R1.update;
R1.cal_posall;

plot_abs_pos(R1.pos_all,alpha,[xA,yA]);


function [r,J]=cal_balance(x,Rod,pA,pdes)
    alpha=pA(3);
    
    alpha_degree=alpha/pi*180;
    
    Rod.theta=x(1:Rod.n_seg);
    Rod.F=x(Rod.n_seg+1:Rod.n_seg+3);
    
    Rod.update;
    
    r=zeros(size(x));
    J=zeros(length(x));
    
    r(1:3)=rotz(alpha_degree)*Rod.pe+pA-pdes;
    r(4:Rod.n_seg+3)=Rod.K_theta*Rod.theta-transpose(Rod.Jacobian)*Rod.F;
    
    J=Rod.cal_rdot;
    J(1:3,1:Rod.n_seg)=rotz(alpha_degree)*Rod.Jacobian;
    J(4:end,1:Rod.n_seg)=Rod.K_theta-Rod.partial;
    J(4:end,end-2:end)=-transpose(Rod.Jacobian);
    
end