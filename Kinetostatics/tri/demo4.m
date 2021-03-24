% 尝试求解顶端约束的Fin-ray型机构，考虑刚性支撑

% 21-03-15 求解效果不好，方程疑似有误，中间的约束力对系统的影响没想清楚。
% 一个也许可行的思路是将约束力平移到末端。然后补上力偶

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


% 参数的选取

% 第一组成功的参数
k1=fix(1/4*nA);
k2=fix(1/4*nB);
Lcon=3/4*sqrt((xA-xB)^2+(yA-yB)^2);
Lcon2=2/4*sqrt((xA-xB)^2+(yA-yB)^2);
LA=0.8*L0;
% 如果使用无刚性约束时的求解结果作为初值，收敛速度会快一点，399→224
% 如果再打开partial项，224→62


% 第二组成功的参数
% k1=fix(3/4*nA);
% k2=fix(3/4*nB);
% Lcon=2/5*sqrt((xA-xB)^2+(yA-yB)^2);
% LA=0.6*L0;
% 惊人的发现：先用无刚性约束时的求解结果作为初值，再打开partial项，收敛速度奇快

% 求解结果很离谱的参数
% k1=fix(1/2*nA);
% k2=fix(1/2*nB);
% Lcon=4/8*sqrt((xA-xB)^2+(yA-yB)^2);
% LA=0.8*L0;
% 关掉partial项，用无刚性约束时的求解结果作为初值，能拿到一个很离谱的图

vA=[ones(k1,1);zeros(nA-k1,1)];
lambdaA=diag(vA);
vB=[ones(k2,1);zeros(nB-k2,1)];
lambdaB=diag(vB);

wid=5e-3;
thi=1e-3;
E=197*1e9;

pdes=[0;0;0];

LB=L0;

RA=planar_nR(E,LA,wid,thi,nA,pdes);
RB=planar_nR(E,LB,wid,thi,nB,pdes);


% 牛顿法求解2N+3方程组

x=zeros(nA+nB+5,1);

%尝试使用无刚性约束时的初值
load('x_init.mat')
x(1:nA+nB+3)=x_init;


TOL=1e-6;
k=1;
while(1)
    if k>500
        error("can not converge")
    end
    
    thetaA=zeros(nA,1);
    thetaB=zeros(nB,1);
    FB=zeros(3,1);
    fcon=0;
    gamma=0;
    
    thetaA=x(1:nA);
    thetaB=x(nA+1:nA+nB);
    FB=x(nA+nB+1:nA+nB+3);
    fcon=x(end-1);
    gamma=x(end);
    
    
    % 开始计算函数值
    
    FA=B*FB;
    
    RA.theta=thetaA;
    RA.F=FA;
    RB.theta=thetaB;
    RB.F=FB;
    
    RA.update;
    RB.update;
    
    pk1=RA.cal_pk(k1);
    pk2=RB.cal_pk(k2);
    
    r2=[pk1(1);pk1(2)]-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*[pk2(1);pk2(2)]+...
        [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[xA-xB+Lcon*cos(gamma);yA-yB+Lcon*sin(gamma)];
    
    r1=RA.pe-A*RB.pe-b;
    r3=RA.K_theta*thetaA-transpose(RA.Jacobian)*FA-lambdaA*transpose(RA.Jacobian)*(-rotz(-alpha_degree))*[fcon*cos(gamma);fcon*sin(gamma);0];
    r4=RB.K_theta*thetaB-transpose(RB.Jacobian)*FB-lambdaB*transpose(RB.Jacobian)*rotz(-beta_degree)*[fcon*cos(gamma);fcon*sin(gamma);0];
    r=[r1;r2;r3;r4];
    
    if(norm(r)<TOL)
        fprintf('Newton Method converge: iteration = %d\n',k-1)
        fprintf('norm(e) = %E\n',norm(r))
        break;
    end
    
    % 开始计算雅可比矩阵
    
    J=zeros(nA+nB+5);
    
    J(1:3,1:nA)=RA.Jacobian;
    J(1:3,nA+1:nA+nB)=-A*RB.Jacobian;
    temp1=RA.Jacobian*lambdaA;
    J(4:5,1:nA)=temp1(1:2,:);
    temp2=RB.Jacobian*lambdaB;
    J(4:5,nA+1:nA+nB)=-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*temp2(1:2,:);
    J(4:5,end)=[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[-Lcon*sin(gamma);Lcon*cos(gamma)];
    
    temp3=RA.partial;
    RA.F=-rotz(-alpha_degree)*[fcon*cos(gamma);fcon*sin(gamma);0];
    RA.cal_partial;
    temp4=lambdaA*RA.partial;
    J(6:nA+5,1:nA)=RA.K_theta-1*(1*temp3+temp4);
    
    J(6:nA+5,nA+nB+1:nA+nB+3)=-transpose(RA.Jacobian)*B;
    J(6:nA+5,end-1)=-lambdaA*transpose(RA.Jacobian)*(-rotz(-alpha_degree)*[cos(gamma);sin(gamma);0]);
    J(6:nA+5,end)=-lambdaA*transpose(RA.Jacobian)*(-rotz(-alpha_degree)*[-fcon*sin(gamma);fcon*cos(gamma);0]);
    
    temp5=RB.partial;
    RB.F=rotz(-beta_degree)*[fcon*cos(gamma);fcon*sin(gamma);0];
    RB.cal_partial;
    temp6=lambdaB*RB.partial;
    J(nA+6:end,nA+1:nA+nB)=RB.K_theta-1*(1*temp5+temp6);
    
    J(nA+6:end,nA+nB+1:nA+nB+3)=-transpose(RB.Jacobian);
    J(nA+6:end,end-1)=-lambdaB*transpose(RB.Jacobian)*(rotz(-beta_degree)*[cos(gamma);sin(gamma);0]);
    J(nA+6:end,end)=-lambdaB*transpose(RB.Jacobian)*(rotz(-beta_degree)*[-fcon*sin(gamma);fcon*cos(gamma);0]);
    
    delta=-pinv(J)*r;
    
    x=x+delta;

    k=k+1;
end

RA.cal_posall;
RB.cal_posall;

A_pos_all=RA.pos_all;
B_pos_all=RB.pos_all;

A_abs_pos_all=plot_abs_pos(A_pos_all,alpha,[xA,yA]);
hold on
B_abs_pos_all=plot_abs_pos(B_pos_all,beta,[xB,yB]);

pk1=RA.cal_pk(k1);
P1(1)=pk1(1)*cos(alpha)-pk1(2)*sin(alpha)+xA;
P1(2)=pk1(1)*sin(alpha)+pk1(2)*cos(alpha)+yA;

pk2=RB.cal_pk(k2);
P2(1)=pk2(1)*cos(beta)-pk2(2)*sin(beta)+xB;
P2(2)=pk2(1)*sin(beta)+pk2(2)*cos(beta)+yB;

plot([P1(1) P2(1)],[P1(2) P2(2)])