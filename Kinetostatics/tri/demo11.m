% 尝试求解顶端约束的Fin-ray型机构
% 考虑内部有4根刚性约束

% 21-03-24晚
% 淦，之前写的是B=-rotz(beta-alpha)，忘了把弧度制转化成角度制了。。。
% 改正过来以后原先不收敛的参数也收敛了...
% 甚至从全为0的初值也能搜索出来了。。。

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
LA=0.8*L0;

% 第1组收敛的参数
% 第1根刚性约束
ka1=fix(1/8*nA);
kb1=fix(1/8*nB);
Lcon1=7/8*sqrt((xA-xB)^2+(yA-yB)^2);

vA1=[ones(ka1,1);zeros(nA-ka1,1)];
lambdaA1=diag(vA1);
vB1=[ones(kb1,1);zeros(nB-kb1,1)];
lambdaB1=diag(vB1);

% 第2根刚性约束
ka2=fix(1/4*nA);
kb2=fix(1/4*nB);
Lcon2=3/4*sqrt((xA-xB)^2+(yA-yB)^2);

vA2=[ones(ka2,1);zeros(nA-ka2,1)];
lambdaA2=diag(vA2);
vB2=[ones(kb2,1);zeros(nB-kb2,1)];
lambdaB2=diag(vB2);

% 第3根刚性约束
ka3=fix(3/8*nA);
kb3=fix(3/8*nB);
Lcon3=5/8*sqrt((xA-xB)^2+(yA-yB)^2);

vA3=[ones(ka3,1);zeros(nA-ka3,1)];
lambdaA3=diag(vA3);
vB3=[ones(kb3,1);zeros(nB-kb3,1)];
lambdaB3=diag(vB3);

% 第4根刚性约束
ka4=fix(4/8*nA);
kb4=fix(4/8*nB);
Lcon4=4/8*sqrt((xA-xB)^2+(yA-yB)^2);

vA4=[ones(ka4,1);zeros(nA-ka4,1)];
lambdaA4=diag(vA4);
vB4=[ones(kb4,1);zeros(nB-kb4,1)];
lambdaB4=diag(vB4);


% 第2组收敛的参数
% % 第1根刚性约束
% ka1=fix(1/5*nA);
% kb1=fix(1/5*nB);
% Lcon1=4/5*sqrt((xA-xB)^2+(yA-yB)^2);
% 
% vA1=[ones(ka1,1);zeros(nA-ka1,1)];
% lambdaA1=diag(vA1);
% vB1=[ones(kb1,1);zeros(nB-kb1,1)];
% lambdaB1=diag(vB1);           
% 
% % 第2根刚性约束
% ka2=fix(2/5*nA);
% kb2=fix(2/5*nB);
% Lcon2=3/5*sqrt((xA-xB)^2+(yA-yB)^2);
% 
% vA2=[ones(ka2,1);zeros(nA-ka2,1)];
% lambdaA2=diag(vA2);
% vB2=[ones(kb2,1);zeros(nB-kb2,1)];
% lambdaB2=diag(vB2);
% 
% % 第3根刚性约束
% ka3=fix(3/5*nA);
% kb3=fix(3/5*nB);
% Lcon3=2/5*sqrt((xA-xB)^2+(yA-yB)^2);
% 
% vA3=[ones(ka3,1);zeros(nA-ka3,1)];
% lambdaA3=diag(vA3);
% vB3=[ones(kb3,1);zeros(nB-kb3,1)];
% lambdaB3=diag(vB3);
% 
% % 第4根刚性约束
% ka4=fix(4/5*nA);
% kb4=fix(4/5*nB);
% Lcon4=1/5*sqrt((xA-xB)^2+(yA-yB)^2);
% 
% vA4=[ones(ka4,1);zeros(nA-ka4,1)];
% lambdaA4=diag(vA4);
% vB4=[ones(kb4,1);zeros(nB-kb4,1)];
% lambdaB4=diag(vB4);


wid=5e-3;
thi=1e-3;
E=197*1e9;

pdes=[0;0;0];

LB=L0;

RA=planar_nR(E,LA,wid,thi,nA,pdes);
RB=planar_nR(E,LB,wid,thi,nB,pdes);


% 牛顿法求解方程组

x=zeros(nA+nB+11,1);

%尝试使用无刚性约束时的初值
% load('x_init.mat')
% x(1:nA+nB+3)=x_init;


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
    fcon3=0;
    gamma3=0;
    fcon4=0;
    gamma4=0;
    
    thetaA=x(1:nA);
    thetaB=x(nA+1:nA+nB);
    FB=x(nA+nB+1:nA+nB+3);
    fcon1=x(nA+nB+4);
    gamma1=x(nA+nB+5);
    fcon2=x(nA+nB+6);
    gamma2=x(nA+nB+7);
    fcon3=x(nA+nB+8);
    gamma3=x(nA+nB+9);
    fcon4=x(nA+nB+10);
    gamma4=x(nA+nB+11);
    
    
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
    
    pka3=RA.cal_pk(ka3);
    pkb3=RB.cal_pk(kb3);
    
    pka4=RA.cal_pk(ka4);
    pkb4=RB.cal_pk(kb4);
    
    r4=[pka1(1);pka1(2)]-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*[pkb1(1);pkb1(2)]+...
        [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[xA-xB+Lcon1*cos(gamma1);yA-yB+Lcon1*sin(gamma1)];
    r5=[pka2(1);pka2(2)]-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*[pkb2(1);pkb2(2)]+...
        [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[xA-xB+Lcon2*cos(gamma2);yA-yB+Lcon2*sin(gamma2)];
    r6=[pka3(1);pka3(2)]-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*[pkb3(1);pkb3(2)]+...
        [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[xA-xB+Lcon3*cos(gamma3);yA-yB+Lcon3*sin(gamma3)];
    r7=[pka4(1);pka4(2)]-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*[pkb4(1);pkb4(2)]+...
        [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[xA-xB+Lcon4*cos(gamma4);yA-yB+Lcon4*sin(gamma4)];
    
    
    r1=RA.pe-A*RB.pe-b;
    r2=RA.K_theta*thetaA-transpose(RA.Jacobian)*FA-...
       lambdaA1*transpose(RA.Jacobian)*(-rotz(-alpha_degree))*[fcon1*cos(gamma1);fcon1*sin(gamma1);0]-...
       lambdaA2*transpose(RA.Jacobian)*(-rotz(-alpha_degree))*[fcon2*cos(gamma2);fcon2*sin(gamma2);0]-...
       lambdaA3*transpose(RA.Jacobian)*(-rotz(-alpha_degree))*[fcon3*cos(gamma3);fcon3*sin(gamma3);0]-...
       lambdaA4*transpose(RA.Jacobian)*(-rotz(-alpha_degree))*[fcon4*cos(gamma4);fcon4*sin(gamma4);0];
    r3=RB.K_theta*thetaB-transpose(RB.Jacobian)*FB-...
       lambdaB1*transpose(RB.Jacobian)*rotz(-beta_degree)*[fcon1*cos(gamma1);fcon1*sin(gamma1);0]-...
       lambdaB2*transpose(RB.Jacobian)*rotz(-beta_degree)*[fcon2*cos(gamma2);fcon2*sin(gamma2);0]-...
       lambdaB3*transpose(RB.Jacobian)*rotz(-beta_degree)*[fcon3*cos(gamma3);fcon3*sin(gamma3);0]-...
       lambdaB4*transpose(RB.Jacobian)*rotz(-beta_degree)*[fcon4*cos(gamma4);fcon4*sin(gamma4);0];
    r=[r1;r2;r3;r4;r5;r6;r7];
    
    norm(r)
    
    if(norm(r)<TOL)
        fprintf('Newton Method converge: iteration = %d\n',k-1)
        fprintf('norm(e) = %E\n',norm(r))
        break;
    end
    
    % 开始计算雅可比矩阵
    
    J=zeros(nA+nB+11);
    
    % 末端约束方程的导数
    J(1:3,1:nA)=RA.Jacobian;
    J(1:3,nA+1:nA+nB)=-A*RB.Jacobian;
    
    
    % 内部约束方程的导数
    % f1
    tempa1=RA.Jacobian*lambdaA1;
    J(nA+nB+4:nA+nB+5,1:nA)=tempa1(1:2,:);
    
    tempb1=RB.Jacobian*lambdaB1;
    J(nA+nB+4:nA+nB+5,nA+1:nA+nB)=-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*tempb1(1:2,:);
    J(nA+nB+4:nA+nB+5,nA+nB+5)=[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[-Lcon1*sin(gamma1);Lcon1*cos(gamma1)];
    
    % f2
    tempa2=RA.Jacobian*lambdaA2;
    J(nA+nB+6:nA+nB+7,1:nA)=tempa2(1:2,:);
    
    tempb2=RB.Jacobian*lambdaB2;
    J(nA+nB+6:nA+nB+7,nA+1:nA+nB)=-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*tempb2(1:2,:);
    J(nA+nB+6:nA+nB+7,nA+nB+7)=[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[-Lcon2*sin(gamma2);Lcon2*cos(gamma2)];
    
    % f3
    tempa3=RA.Jacobian*lambdaA3;
    J(nA+nB+8:nA+nB+9,1:nA)=tempa3(1:2,:);
    
    tempb3=RB.Jacobian*lambdaB3;
    J(nA+nB+8:nA+nB+9,nA+1:nA+nB)=-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*tempb3(1:2,:);
    J(nA+nB+8:nA+nB+9,nA+nB+9)=[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[-Lcon3*sin(gamma3);Lcon3*cos(gamma3)];
    
    % f4
    tempa4=RA.Jacobian*lambdaA4;
    J(nA+nB+10:nA+nB+11,1:nA)=tempa4(1:2,:);
    
    tempb4=RB.Jacobian*lambdaB4;
    J(nA+nB+10:nA+nB+11,nA+1:nA+nB)=-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*tempb4(1:2,:);
    J(nA+nB+10:nA+nB+11,nA+nB+11)=[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[-Lcon4*sin(gamma4);Lcon4*cos(gamma4)];
    
    
    % 受力平衡方程的导数
    % fa
    % fa对theta_a的导数
    % 求和的第0项 对应Fb
    tempA0=RA.partial;
    % 求和的第1项 对应fcon1
    RA.F=-rotz(-alpha_degree)*[fcon1*cos(gamma1);fcon1*sin(gamma1);0];
    RA.cal_partial;
    tempA1=lambdaA1*RA.partial;
    % 求和的第2项 对应fcon2
    RA.F=-rotz(-alpha_degree)*[fcon2*cos(gamma2);fcon2*sin(gamma2);0];
    RA.cal_partial;
    tempA2=lambdaA2*RA.partial;
    % 求和的第3项 对应fcon3
    RA.F=-rotz(-alpha_degree)*[fcon3*cos(gamma3);fcon3*sin(gamma3);0];
    RA.cal_partial;
    tempA3=lambdaA3*RA.partial;
    % 求和的第4项 对应fcon4
    RA.F=-rotz(-alpha_degree)*[fcon4*cos(gamma4);fcon4*sin(gamma4);0];
    RA.cal_partial;
    tempA4=lambdaA4*RA.partial;
    
    J(4:nA+3,1:nA)=RA.K_theta-1*(1*tempA0+tempA1+tempA2+tempA3+tempA4);
    
    % fa对F_b的导数
    J(4:nA+3,nA+nB+1:nA+nB+3)=-transpose(RA.Jacobian)*B;
    
    % fa对第1组约束的导数
    J(4:nA+3,nA+nB+4)=-lambdaA1*transpose(RA.Jacobian)*(-rotz(-alpha_degree)*[cos(gamma1);sin(gamma1);0]);
    J(4:nA+3,nA+nB+5)=-lambdaA1*transpose(RA.Jacobian)*(-rotz(-alpha_degree)*[-fcon1*sin(gamma1);fcon1*cos(gamma1);0]);
    % fa对第2组约束的导数
    J(4:nA+3,nA+nB+6)=-lambdaA2*transpose(RA.Jacobian)*(-rotz(-alpha_degree)*[cos(gamma2);sin(gamma2);0]);
    J(4:nA+3,nA+nB+7)=-lambdaA2*transpose(RA.Jacobian)*(-rotz(-alpha_degree)*[-fcon2*sin(gamma2);fcon1*cos(gamma2);0]);
    % fa对第3组约束的导数
    J(4:nA+3,nA+nB+8)=-lambdaA3*transpose(RA.Jacobian)*(-rotz(-alpha_degree)*[cos(gamma3);sin(gamma3);0]);
    J(4:nA+3,nA+nB+9)=-lambdaA3*transpose(RA.Jacobian)*(-rotz(-alpha_degree)*[-fcon3*sin(gamma3);fcon3*cos(gamma3);0]);
    % fa对第4组约束的导数
    J(4:nA+3,nA+nB+10)=-lambdaA4*transpose(RA.Jacobian)*(-rotz(-alpha_degree)*[cos(gamma4);sin(gamma4);0]);
    J(4:nA+3,nA+nB+11)=-lambdaA4*transpose(RA.Jacobian)*(-rotz(-alpha_degree)*[-fcon4*sin(gamma4);fcon4*cos(gamma4);0]);
    
    % fb
    % fb对theta_b的导数
    % 求和的第0项 对应Fb
    tempB0=RB.partial;
    % 求和的第1项 对应fcon1
    RB.F=rotz(-beta_degree)*[fcon1*cos(gamma1);fcon1*sin(gamma1);0];
    RB.cal_partial;
    tempB1=lambdaB1*RB.partial;
    % 求和的第2项 对应fcon2
    RB.F=rotz(-beta_degree)*[fcon2*cos(gamma2);fcon2*sin(gamma2);0];
    RB.cal_partial;
    tempB2=lambdaB2*RB.partial;
    % 求和的第3项 对应fcon3
    RB.F=rotz(-beta_degree)*[fcon3*cos(gamma3);fcon3*sin(gamma3);0];
    RB.cal_partial;
    tempB3=lambdaB3*RB.partial;
    % 求和的第4项 对应fcon4
    RB.F=rotz(-beta_degree)*[fcon4*cos(gamma4);fcon4*sin(gamma4);0];
    RB.cal_partial;
    tempB4=lambdaB4*RB.partial;
    
    J(nA+4:nA+nB+3,nA+1:nA+nB)=RB.K_theta-1*(1*tempB0+tempB1+tempB2+tempB3+tempB4);
    
    % fb对F_b的导数
    J(nA+4:nA+nB+3,nA+nB+1:nA+nB+3)=-transpose(RB.Jacobian);
    
    % fb对第1组约束的导数
    J(nA+4:nA+nB+3,nA+nB+4)=-lambdaB1*transpose(RB.Jacobian)*(rotz(-beta_degree)*[cos(gamma1);sin(gamma1);0]);
    J(nA+4:nA+nB+3,nA+nB+5)=-lambdaB1*transpose(RB.Jacobian)*(rotz(-beta_degree)*[-fcon1*sin(gamma1);fcon1*cos(gamma1);0]);
    % fb对第2组约束的导数
    J(nA+4:nA+nB+3,nA+nB+6)=-lambdaB2*transpose(RB.Jacobian)*(rotz(-beta_degree)*[cos(gamma2);sin(gamma2);0]);
    J(nA+4:nA+nB+3,nA+nB+7)=-lambdaB2*transpose(RB.Jacobian)*(rotz(-beta_degree)*[-fcon2*sin(gamma2);fcon2*cos(gamma2);0]);
    % fb对第3组约束的导数
    J(nA+4:nA+nB+3,nA+nB+8)=-lambdaB3*transpose(RB.Jacobian)*(rotz(-beta_degree)*[cos(gamma3);sin(gamma3);0]);
    J(nA+4:nA+nB+3,nA+nB+9)=-lambdaB3*transpose(RB.Jacobian)*(rotz(-beta_degree)*[-fcon3*sin(gamma3);fcon3*cos(gamma3);0]);
    % fb对第4组约束的导数
    J(nA+4:nA+nB+3,nA+nB+10)=-lambdaB4*transpose(RB.Jacobian)*(rotz(-beta_degree)*[cos(gamma4);sin(gamma4);0]);
    J(nA+4:nA+nB+3,nA+nB+11)=-lambdaB4*transpose(RB.Jacobian)*(rotz(-beta_degree)*[-fcon4*sin(gamma4);fcon4*cos(gamma4);0]);
    
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


% 画出内部的刚性约束
% 第1根刚性约束
pka1=RA.cal_pk(ka1);
PA1(1)=pka1(1)*cos(alpha)-pka1(2)*sin(alpha)+xA;
PA1(2)=pka1(1)*sin(alpha)+pka1(2)*cos(alpha)+yA;

pkb1=RB.cal_pk(kb1);
PB1(1)=pkb1(1)*cos(beta)-pkb1(2)*sin(beta)+xB;
PB1(2)=pkb1(1)*sin(beta)+pkb1(2)*cos(beta)+yB;

plot([PA1(1) PB1(1)],[PA1(2) PB1(2)])

% 第2根刚性约束
pka2=RA.cal_pk(ka2);
PA2(1)=pka2(1)*cos(alpha)-pka2(2)*sin(alpha)+xA;
PA2(2)=pka2(1)*sin(alpha)+pka2(2)*cos(alpha)+yA;

pkb2=RB.cal_pk(kb2);
PB2(1)=pkb2(1)*cos(beta)-pkb2(2)*sin(beta)+xB;
PB2(2)=pkb2(1)*sin(beta)+pkb2(2)*cos(beta)+yB;

plot([PA2(1) PB2(1)],[PA2(2) PB2(2)])

% 第3根刚性约束
pka3=RA.cal_pk(ka3);
PA3(1)=pka3(1)*cos(alpha)-pka3(2)*sin(alpha)+xA;
PA3(2)=pka3(1)*sin(alpha)+pka3(2)*cos(alpha)+yA;

pkb3=RB.cal_pk(kb3);
PB3(1)=pkb3(1)*cos(beta)-pkb3(2)*sin(beta)+xB;
PB3(2)=pkb3(1)*sin(beta)+pkb3(2)*cos(beta)+yB;

plot([PA3(1) PB3(1)],[PA3(2) PB3(2)])

% 第4根刚性约束
pka4=RA.cal_pk(ka4);
PA4(1)=pka4(1)*cos(alpha)-pka4(2)*sin(alpha)+xA;
PA4(2)=pka4(1)*sin(alpha)+pka4(2)*cos(alpha)+yA;

pkb4=RB.cal_pk(kb4);
PB4(1)=pkb4(1)*cos(beta)-pkb4(2)*sin(beta)+xB;
PB4(2)=pkb4(1)*sin(beta)+pkb4(2)*cos(beta)+yB;

plot([PA4(1) PB4(1)],[PA4(2) PB4(2)])