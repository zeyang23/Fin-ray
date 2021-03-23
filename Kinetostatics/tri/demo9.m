% 尝试求解顶端约束的Fin-ray型机构
% 考虑内部有两根刚性约束

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

B=-rotz(beta-alpha);


nA=50;
nB=50;


% 参数的选取
LA=0.8*L0;

% 第1根刚性约束
ka1=fix(1/4*nA);
kb1=fix(1/4*nB);
Lcon1=3/4*sqrt((xA-xB)^2+(yA-yB)^2);

vA1=[ones(ka1,1);zeros(nA-ka1,1)];
lambdaA1=diag(vA1);
vB1=[ones(kb1,1);zeros(nB-kb1,1)];
lambdaB1=diag(vB1);


% 第2根刚性约束

% 第1组收敛参数
% ka2=fix(3/4*nA);
% kb2=fix(3/4*nB);
% Lcon2=1/4*sqrt((xA-xB)^2+(yA-yB)^2);

% 第2组收敛参数
ka2=fix(1/8*nA);
kb2=fix(1/8*nB);
Lcon2=7/8*sqrt((xA-xB)^2+(yA-yB)^2);

vA2=[ones(ka2,1);zeros(nA-ka2,1)];
lambdaA2=diag(vA2);
vB2=[ones(kb2,1);zeros(nB-kb2,1)];
lambdaB2=diag(vB2);



wid=5e-3;
thi=1e-3;
E=197*1e9;

pdes=[0;0;0];

LB=L0;

RA=planar_nR(E,LA,wid,thi,nA,pdes);
RB=planar_nR(E,LB,wid,thi,nB,pdes);


% 牛顿法求解方程组

x=zeros(nA+nB+7,1);

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
    fcon1=0;
    gamma1=0;
    fcon2=0;
    gamma2=0;
    
    thetaA=x(1:nA);
    thetaB=x(nA+1:nA+nB);
    FB=x(nA+nB+1:nA+nB+3);
    fcon1=x(end-3);
    gamma1=x(end-2);
    fcon2=x(end-1);
    gamma2=x(end);
    
    
    % 开始计算函数值
    
    FA=B*FB;
    
    RA.theta=thetaA;
    RA.F=FA;
    RB.theta=thetaB;
    RB.F=FB;
    
    RA.update;
    RB.update;
    
    pka1=RA.cal_pk(ka1);
    pkb1=RB.cal_pk(kb1);
    
    pka2=RA.cal_pk(ka2);
    pkb2=RB.cal_pk(kb2);
    
    r4=[pka1(1);pka1(2)]-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*[pkb1(1);pkb1(2)]+...
        [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[xA-xB+Lcon1*cos(gamma1);yA-yB+Lcon1*sin(gamma1)];
    r5=[pka2(1);pka2(2)]-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*[pkb2(1);pkb2(2)]+...
        [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[xA-xB+Lcon2*cos(gamma2);yA-yB+Lcon2*sin(gamma2)];
    
    
    r1=RA.pe-A*RB.pe-b;
    r2=RA.K_theta*thetaA-transpose(RA.Jacobian)*FA-lambdaA1*transpose(RA.Jacobian)*(-rotz(-alpha_degree))*[fcon1*cos(gamma1);fcon1*sin(gamma1);0]-...
       lambdaA2*transpose(RA.Jacobian)*(-rotz(-alpha_degree))*[fcon2*cos(gamma2);fcon2*sin(gamma2);0];
    r3=RB.K_theta*thetaB-transpose(RB.Jacobian)*FB-lambdaB1*transpose(RB.Jacobian)*rotz(-beta_degree)*[fcon1*cos(gamma1);fcon1*sin(gamma1);0]-...
       lambdaB2*transpose(RB.Jacobian)*rotz(-beta_degree)*[fcon2*cos(gamma2);fcon2*sin(gamma2);0];
    r=[r1;r2;r3;r4;r5];
    
    norm(r)
    
    if(norm(r)<TOL)
        fprintf('Newton Method converge: iteration = %d\n',k-1)
        fprintf('norm(e) = %E\n',norm(r))
        break;
    end
    
    % 开始计算雅可比矩阵
    
    J=zeros(nA+nB+7);
    
    J(1:3,1:nA)=RA.Jacobian;
    J(1:3,nA+1:nA+nB)=-A*RB.Jacobian;
    
    
    tempa1=RA.Jacobian*lambdaA1;
    J(end-3:end-2,1:nA)=tempa1(1:2,:);
    
    tempb1=RB.Jacobian*lambdaB1;
    J(end-3:end-2,nA+1:nA+nB)=-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*tempb1(1:2,:);
    J(end-3:end-2,end-2)=[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[-Lcon1*sin(gamma1);Lcon1*cos(gamma1)];
    
    tempa2=RA.Jacobian*lambdaA2;
    J(end-1:end,1:nA)=tempa2(1:2,:);
    
    tempb2=RB.Jacobian*lambdaB2;
    J(end-1:end,nA+1:nA+nB)=-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*tempb2(1:2,:);
    J(end-1:end,end)=[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[-Lcon2*sin(gamma2);Lcon2*cos(gamma2)];
    
    
    temp3=RA.partial;
    RA.F=-rotz(-alpha_degree)*[fcon1*cos(gamma1);fcon1*sin(gamma1);0];
    RA.cal_partial;
    temp4=lambdaA1*RA.partial;
    RA.F=-rotz(-alpha_degree)*[fcon2*cos(gamma2);fcon2*sin(gamma2);0];
    RA.cal_partial;
    temp5=lambdaA2*RA.partial;
    
    J(4:nA+3,1:nA)=RA.K_theta-1*(1*temp3+temp4+temp5);
    
    J(4:nA+3,nA+nB+1:nA+nB+3)=-transpose(RA.Jacobian)*B;
    J(4:nA+3,end-3)=-lambdaA1*transpose(RA.Jacobian)*(-rotz(-alpha_degree)*[cos(gamma1);sin(gamma1);0]);
    J(4:nA+3,end-2)=-lambdaA1*transpose(RA.Jacobian)*(-rotz(-alpha_degree)*[-fcon1*sin(gamma1);fcon1*cos(gamma1);0]);
    J(4:nA+3,end-1)=-lambdaA2*transpose(RA.Jacobian)*(-rotz(-alpha_degree)*[cos(gamma2);sin(gamma2);0]);
    J(4:nA+3,end)=-lambdaA2*transpose(RA.Jacobian)*(-rotz(-alpha_degree)*[-fcon2*sin(gamma2);fcon1*cos(gamma2);0]);
    
    
    temp6=RB.partial;
    RB.F=rotz(-beta_degree)*[fcon1*cos(gamma1);fcon1*sin(gamma1);0];
    RB.cal_partial;
    temp7=lambdaB1*RB.partial;
    RB.F=rotz(-beta_degree)*[fcon2*cos(gamma2);fcon2*sin(gamma2);0];
    RB.cal_partial;
    temp8=lambdaB2*RB.partial;
    
    J(nA+4:nA+nB+3,nA+1:nA+nB)=RB.K_theta-1*(1*temp6+temp7+temp8);
    
    J(nA+4:nA+nB+3,nA+nB+1:nA+nB+3)=-transpose(RB.Jacobian);
    J(nA+4:nA+nB+3,end-3)=-lambdaB1*transpose(RB.Jacobian)*(rotz(-beta_degree)*[cos(gamma1);sin(gamma1);0]);
    J(nA+4:nA+nB+3,end-2)=-lambdaB1*transpose(RB.Jacobian)*(rotz(-beta_degree)*[-fcon1*sin(gamma1);fcon1*cos(gamma1);0]);
    J(nA+4:nA+nB+3,end-1)=-lambdaB2*transpose(RB.Jacobian)*(rotz(-beta_degree)*[cos(gamma2);sin(gamma2);0]);
    J(nA+4:nA+nB+3,end)=-lambdaB2*transpose(RB.Jacobian)*(rotz(-beta_degree)*[-fcon2*sin(gamma2);fcon2*cos(gamma2);0]);
    
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

pka1=RA.cal_pk(ka1);
P1(1)=pka1(1)*cos(alpha)-pka1(2)*sin(alpha)+xA;
P1(2)=pka1(1)*sin(alpha)+pka1(2)*cos(alpha)+yA;

pkb1=RB.cal_pk(kb1);
P2(1)=pkb1(1)*cos(beta)-pkb1(2)*sin(beta)+xB;
P2(2)=pkb1(1)*sin(beta)+pkb1(2)*cos(beta)+yB;

plot([P1(1) P2(1)],[P1(2) P2(2)])

pka2=RA.cal_pk(ka2);
P3(1)=pka2(1)*cos(alpha)-pka2(2)*sin(alpha)+xA;
P3(2)=pka2(1)*sin(alpha)+pka2(2)*cos(alpha)+yA;

pkb2=RB.cal_pk(kb2);
P4(1)=pkb2(1)*cos(beta)-pkb2(2)*sin(beta)+xB;
P4(2)=pkb2(1)*sin(beta)+pkb2(2)*cos(beta)+yB;

plot([P3(1) P4(1)],[P3(2) P4(2)])