% 求解顶端约束的Fin-ray型机构，无刚性支撑
% 使用fsolve

% GradientCheck能够通过。注意partial项要打开

% 21-03-24晚
% 淦，之前写的是B=-rotz(beta-alpha)，忘了把弧度制转化成角度制了。。。

clear
clc

L0=1;

xA=0;
yA=0;

xB=0.35*L0;
yB=0*L0;

psi_degree=20;
alpha_degree=80;
beta_degree=100;

psi=psi_degree/180*pi;
alpha=alpha_degree/180*pi;
beta=beta_degree/180*pi;

A=rotz(-alpha_degree)*rotz(beta_degree);
b=rotz(-alpha_degree)*[xB-xA;yB-yA;beta-alpha-psi];

B=-rotz(beta_degree-alpha_degree);


nA=50;
nB=50;

wid=5e-3;
thi=1e-3;
E=197*1e9;

pdes=[0;0;0];

LA=0.8*L0;
LB=L0;

RA=planar_nR(E,LA,wid,thi,nA,pdes);
RB=planar_nR(E,LB,wid,thi,nB,pdes);


% 牛顿法求解2N+3方程组
x0=zeros(nA+nB+3,1);

options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',true);

f=@(x) myfun(x,nA,nB,A,B,b,RA,RB);
[x,fval,exitflag,output] = fsolve(f,x0,options);

RA.cal_posall;
RB.cal_posall;

A_pos_all=RA.pos_all;
B_pos_all=RB.pos_all;

A_abs_pos_all=plot_abs_pos(A_pos_all,alpha,[xA,yA]);
hold on
B_abs_pos_all=plot_abs_pos(B_pos_all,beta,[xB,yB]);


function [r,J]=myfun(x,nA,nB,A,B,b,RA,RB)
    
    thetaA=zeros(nA,1);
    thetaB=zeros(nB,1);
    FB=zeros(3,1);
    
    thetaA=x(1:nA);
    thetaB=x(nA+1:nA+nB);
    FB=x(end-2:end);
    
    FA=B*FB;
    
    RA.theta=thetaA;
    RA.F=FA;
    RB.theta=thetaB;
    RB.F=FB;
    
    RA.update;
    RB.update;
    
    r1=RA.pe-A*RB.pe-b;
    r2=RA.K_theta*thetaA-transpose(RA.Jacobian)*FA;
    r3=RB.K_theta*thetaB-transpose(RB.Jacobian)*FB;
    r=[r1;r2;r3];
    
    if nargout > 1
        J=zeros(nA+nB+3);

        J(1:3,1:nA)=RA.Jacobian;
        J(1:3,nA+1:nA+nB)=-A*RB.Jacobian;
        J(4:nA+3,1:nA)=RA.K_theta-1*RA.partial;
        J(4:nA+3,end-2:end)=-transpose(RA.Jacobian)*B;
        J(nA+4:end,nA+1:nA+nB)=RB.K_theta-1*RB.partial;
        J(nA+4:end,end-2:end)=-transpose(RB.Jacobian);
    end
end