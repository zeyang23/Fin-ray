% 使用fsolve求解
% 考虑内部有4根刚性约束

% 21-03-24
% 无法通过GradientCheck 梯度的计算可能有问题
% 不给fsolve提供梯度也能解出来，只不过速度很慢

% fsolve好强。。。

% 21-03-24晚
% 淦，之前写的是B=-rotz(beta-alpha)，忘了把弧度制转化成角度制了。。。

% 21-03-26
% 终于能把梯度算对了。可以通过GradientCheck
% 值得注意的是，('FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-10)这组检查参数是无法通过的。
% 事实上，之前的无刚性约束的方程在这组检查参数下也无法通过。这个事情也许值得思考。

% 在梯度终于算对了以后，算法的迭代次数显著减小，鲁棒性显著增强


clear
clc

%% 参数设置
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
% % 第1根刚性约束
% ka1=fix(1/8*nA);
% kb1=fix(1/8*nB);
% Lcon1=7/8*sqrt((xA-xB)^2+(yA-yB)^2);
% 
% vA1=[ones(ka1,1);zeros(nA-ka1,1)];
% lambdaA1=diag(vA1);
% vB1=[ones(kb1,1);zeros(nB-kb1,1)];
% lambdaB1=diag(vB1);
% 
% % 第2根刚性约束
% ka2=fix(1/4*nA);
% kb2=fix(1/4*nB);
% Lcon2=3/4*sqrt((xA-xB)^2+(yA-yB)^2);
% 
% vA2=[ones(ka2,1);zeros(nA-ka2,1)];
% lambdaA2=diag(vA2);
% vB2=[ones(kb2,1);zeros(nB-kb2,1)];
% lambdaB2=diag(vB2);
% 
% % 第3根刚性约束
% ka3=fix(3/8*nA);
% kb3=fix(3/8*nB);
% Lcon3=5/8*sqrt((xA-xB)^2+(yA-yB)^2);
% 
% vA3=[ones(ka3,1);zeros(nA-ka3,1)];
% lambdaA3=diag(vA3);
% vB3=[ones(kb3,1);zeros(nB-kb3,1)];
% lambdaB3=diag(vB3);
% 
% % 第4根刚性约束
% ka4=fix(4/8*nA);
% kb4=fix(4/8*nB);
% Lcon4=4/8*sqrt((xA-xB)^2+(yA-yB)^2);
% 
% vA4=[ones(ka4,1);zeros(nA-ka4,1)];
% lambdaA4=diag(vA4);
% vB4=[ones(kb4,1);zeros(nB-kb4,1)];
% lambdaB4=diag(vB4);


% 第2组收敛的参数
% 第1根刚性约束
ka1=fix(1/5*nA);
kb1=fix(1/5*nB);
Lcon1=4/5*sqrt((xA-xB)^2+(yA-yB)^2);

LambdaA1=zeros(ka1,nA);
LambdaA1(1:ka1,1:ka1)=diag(ones(ka1,1));
LambdaB1=zeros(kb1,nB);
LambdaB1(1:kb1,1:kb1)=diag(ones(kb1,1));          

% 第2根刚性约束
ka2=fix(2/5*nA);
kb2=fix(2/5*nB);
Lcon2=3/5*sqrt((xA-xB)^2+(yA-yB)^2);

LambdaA2=zeros(ka2,nA);
LambdaA2(1:ka2,1:ka2)=diag(ones(ka2,1));
LambdaB2=zeros(kb2,nB);
LambdaB2(1:kb2,1:kb2)=diag(ones(kb2,1));   

% 第3根刚性约束
ka3=fix(3/5*nA);
kb3=fix(3/5*nB);
Lcon3=2/5*sqrt((xA-xB)^2+(yA-yB)^2);

LambdaA3=zeros(ka3,nA);
LambdaA3(1:ka3,1:ka3)=diag(ones(ka3,1));
LambdaB3=zeros(kb3,nB);
LambdaB3(1:kb3,1:kb3)=diag(ones(kb3,1));   

% 第4根刚性约束
ka4=fix(4/5*nA);
kb4=fix(4/5*nB);
Lcon4=1/5*sqrt((xA-xB)^2+(yA-yB)^2);

LambdaA4=zeros(ka4,nA);
LambdaA4(1:ka4,1:ka4)=diag(ones(ka4,1));
LambdaB4=zeros(kb4,nB);
LambdaB4(1:kb4,1:kb4)=diag(ones(kb4,1));   


wid=5e-3;
thi=1e-3;
E=197*1e9;

pdes=[0;0;0];

LB=L0;

RA=planar_nR(E,LA,wid,thi,nA,pdes);
RAk1=planar_nR(E,ka1/nA*LA,wid,thi,ka1,pdes);
RAk2=planar_nR(E,ka2/nA*LA,wid,thi,ka2,pdes);
RAk3=planar_nR(E,ka3/nA*LA,wid,thi,ka3,pdes);
RAk4=planar_nR(E,ka4/nA*LA,wid,thi,ka4,pdes);

RB=planar_nR(E,LB,wid,thi,nB,pdes);
RBk1=planar_nR(E,kb1/nB*LB,wid,thi,kb1,pdes);
RBk2=planar_nR(E,kb2/nB*LB,wid,thi,kb2,pdes);
RBk3=planar_nR(E,kb3/nB*LB,wid,thi,kb3,pdes);
RBk4=planar_nR(E,kb4/nB*LB,wid,thi,kb4,pdes);


%% fsolve求解方程组

x0=zeros(nA+nB+11,1);

%尝试使用无刚性约束时的初值
% load('x_init.mat')
% x0(1:nA+nB+3)=x_init;

options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false);
% options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',true,'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-10);

f=@(x) myfun(x,xA,xB,yA,yB,beta,beta_degree,alpha,alpha_degree,A,B,b,nA,nB,RA,RB,RAk1,RAk2,RAk3,RAk4,RBk1,RBk2,RBk3,RBk4,Lcon1,Lcon2,Lcon3,Lcon4,ka1,ka2,ka3,ka4,kb1,kb2,kb3,kb4,LambdaA1,LambdaA2,LambdaA3,LambdaA4,LambdaB1,LambdaB2,LambdaB3,LambdaB4);
[x,fval,exitflag,output] = fsolve(f,x0,options);


%% 画出内部的刚性约束
RA.cal_posall;
RB.cal_posall;

A_pos_all=RA.pos_all;
B_pos_all=RB.pos_all;

A_abs_pos_all=plot_abs_pos(A_pos_all,alpha,[xA,yA]);
hold on
B_abs_pos_all=plot_abs_pos(B_pos_all,beta,[xB,yB]);


% 第1根刚性约束
pka1=RAk1.pe;
PA1(1)=pka1(1)*cos(alpha)-pka1(2)*sin(alpha)+xA;
PA1(2)=pka1(1)*sin(alpha)+pka1(2)*cos(alpha)+yA;

pkb1=RBk1.pe;
PB1(1)=pkb1(1)*cos(beta)-pkb1(2)*sin(beta)+xB;
PB1(2)=pkb1(1)*sin(beta)+pkb1(2)*cos(beta)+yB;

plot([PA1(1) PB1(1)],[PA1(2) PB1(2)])

% 第2根刚性约束
pka2=RAk2.pe;
PA2(1)=pka2(1)*cos(alpha)-pka2(2)*sin(alpha)+xA;
PA2(2)=pka2(1)*sin(alpha)+pka2(2)*cos(alpha)+yA;

pkb2=RBk2.pe;
PB2(1)=pkb2(1)*cos(beta)-pkb2(2)*sin(beta)+xB;
PB2(2)=pkb2(1)*sin(beta)+pkb2(2)*cos(beta)+yB;

plot([PA2(1) PB2(1)],[PA2(2) PB2(2)])

% 第3根刚性约束
pka3=RAk3.pe;
PA3(1)=pka3(1)*cos(alpha)-pka3(2)*sin(alpha)+xA;
PA3(2)=pka3(1)*sin(alpha)+pka3(2)*cos(alpha)+yA;

pkb3=RBk3.pe;
PB3(1)=pkb3(1)*cos(beta)-pkb3(2)*sin(beta)+xB;
PB3(2)=pkb3(1)*sin(beta)+pkb3(2)*cos(beta)+yB;

plot([PA3(1) PB3(1)],[PA3(2) PB3(2)])

% 第4根刚性约束
pka4=RAk4.pe;
PA4(1)=pka4(1)*cos(alpha)-pka4(2)*sin(alpha)+xA;
PA4(2)=pka4(1)*sin(alpha)+pka4(2)*cos(alpha)+yA;

pkb4=RBk4.pe;
PB4(1)=pkb4(1)*cos(beta)-pkb4(2)*sin(beta)+xB;
PB4(2)=pkb4(1)*sin(beta)+pkb4(2)*cos(beta)+yB;

plot([PA4(1) PB4(1)],[PA4(2) PB4(2)])


%%
function [r,J]=myfun(x,xA,xB,yA,yB,beta,beta_degree,alpha,alpha_degree,A,B,b,nA,nB,RA,RB,RAk1,RAk2,RAk3,RAk4,RBk1,RBk2,RBk3,RBk4,Lcon1,Lcon2,Lcon3,Lcon4,ka1,ka2,ka3,ka4,kb1,kb2,kb3,kb4,lambdaA1,lambdaA2,lambdaA3,lambdaA4,lambdaB1,lambdaB2,lambdaB3,lambdaB4)
    
    thetaA=zeros(nA,1);
    thetaAk1=zeros(ka1,1);
    thetaAk2=zeros(ka2,1);
    thetaAk3=zeros(ka3,1);
    thetaAk4=zeros(ka4,1);
    
    thetaB=zeros(nB,1);
    thetaBk1=zeros(kb1,1);
    thetaBk2=zeros(kb2,1);
    thetaBk3=zeros(kb3,1);
    thetaBk4=zeros(kb4,1);
    
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
    thetaAk1=thetaA(1:ka1);
    thetaAk2=thetaA(1:ka2);
    thetaAk3=thetaA(1:ka3);
    thetaAk4=thetaA(1:ka4);
    
    thetaB=x(nA+1:nA+nB);
    thetaBk1=thetaB(1:kb1);
    thetaBk2=thetaB(1:kb2);
    thetaBk3=thetaB(1:kb3);
    thetaBk4=thetaB(1:kb4);
    
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
    
    RAk1.theta=thetaAk1;
    RAk1.F=(-rotz(-alpha_degree))*[fcon1*cos(gamma1);fcon1*sin(gamma1);0];
    RAk2.theta=thetaAk2;
    RAk2.F=(-rotz(-alpha_degree))*[fcon2*cos(gamma2);fcon2*sin(gamma2);0];
    RAk3.theta=thetaAk3;
    RAk3.F=(-rotz(-alpha_degree))*[fcon3*cos(gamma3);fcon3*sin(gamma3);0];
    RAk4.theta=thetaAk4;
    RAk4.F=(-rotz(-alpha_degree))*[fcon4*cos(gamma4);fcon4*sin(gamma4);0];
    
    RBk1.theta=thetaBk1;
    RBk1.F=rotz(-beta_degree)*[fcon1*cos(gamma1);fcon1*sin(gamma1);0];
    RBk2.theta=thetaBk2;
    RBk2.F=rotz(-beta_degree)*[fcon2*cos(gamma2);fcon2*sin(gamma2);0];
    RBk3.theta=thetaBk3;
    RBk3.F=rotz(-beta_degree)*[fcon3*cos(gamma3);fcon3*sin(gamma3);0];
    RBk4.theta=thetaBk4;
    RBk4.F=rotz(-beta_degree)*[fcon4*cos(gamma4);fcon4*sin(gamma4);0];
    
    
    RA.update;
    RAk1.update;
    RAk2.update;
    RAk3.update;
    RAk4.update;
    
    RB.update;
    RBk1.update;
    RBk2.update;
    RBk3.update;
    RBk4.update;
    
    
    pka1=RAk1.pe;
    pkb1=RBk1.pe;
    
    pka2=RAk2.pe;
    pkb2=RBk2.pe;
    
    pka3=RAk3.pe;
    pkb3=RBk3.pe;
    
    pka4=RAk4.pe;
    pkb4=RBk4.pe;
    
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
       transpose(lambdaA1)*transpose(RAk1.Jacobian)*(-rotz(-alpha_degree))*[fcon1*cos(gamma1);fcon1*sin(gamma1);0]-...
       transpose(lambdaA2)*transpose(RAk2.Jacobian)*(-rotz(-alpha_degree))*[fcon2*cos(gamma2);fcon2*sin(gamma2);0]-...
       transpose(lambdaA3)*transpose(RAk3.Jacobian)*(-rotz(-alpha_degree))*[fcon3*cos(gamma3);fcon3*sin(gamma3);0]-...
       transpose(lambdaA4)*transpose(RAk4.Jacobian)*(-rotz(-alpha_degree))*[fcon4*cos(gamma4);fcon4*sin(gamma4);0];
    r3=RB.K_theta*thetaB-transpose(RB.Jacobian)*FB-...
       transpose(lambdaB1)*transpose(RBk1.Jacobian)*rotz(-beta_degree)*[fcon1*cos(gamma1);fcon1*sin(gamma1);0]-...
       transpose(lambdaB2)*transpose(RBk2.Jacobian)*rotz(-beta_degree)*[fcon2*cos(gamma2);fcon2*sin(gamma2);0]-...
       transpose(lambdaB3)*transpose(RBk3.Jacobian)*rotz(-beta_degree)*[fcon3*cos(gamma3);fcon3*sin(gamma3);0]-...
       transpose(lambdaB4)*transpose(RBk4.Jacobian)*rotz(-beta_degree)*[fcon4*cos(gamma4);fcon4*sin(gamma4);0];
    r=[r1;r2;r3;r4;r5;r6;r7];

    
    % 开始计算雅可比矩阵
    if nargout > 1
        J=zeros(nA+nB+11);

        % 末端约束方程的导数
        J(1:3,1:nA)=RA.Jacobian;
        J(1:3,nA+1:nA+nB)=-A*RB.Jacobian;


        % 内部约束方程的导数
        % f1
        tempa1=RAk1.Jacobian*lambdaA1;
        J(nA+nB+4:nA+nB+5,1:nA)=tempa1(1:2,:);

        tempb1=RBk1.Jacobian*lambdaB1;
        J(nA+nB+4:nA+nB+5,nA+1:nA+nB)=-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*tempb1(1:2,:);
        J(nA+nB+4:nA+nB+5,nA+nB+5)=[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[-Lcon1*sin(gamma1);Lcon1*cos(gamma1)];

        % f2
        tempa2=RAk2.Jacobian*lambdaA2;
        J(nA+nB+6:nA+nB+7,1:nA)=tempa2(1:2,:);

        tempb2=RBk2.Jacobian*lambdaB2;
        J(nA+nB+6:nA+nB+7,nA+1:nA+nB)=-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*tempb2(1:2,:);
        J(nA+nB+6:nA+nB+7,nA+nB+7)=[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[-Lcon2*sin(gamma2);Lcon2*cos(gamma2)];

        % f3
        tempa3=RAk3.Jacobian*lambdaA3;
        J(nA+nB+8:nA+nB+9,1:nA)=tempa3(1:2,:);

        tempb3=RBk3.Jacobian*lambdaB3;
        J(nA+nB+8:nA+nB+9,nA+1:nA+nB)=-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*tempb3(1:2,:);
        J(nA+nB+8:nA+nB+9,nA+nB+9)=[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[-Lcon3*sin(gamma3);Lcon3*cos(gamma3)];

        % f4
        tempa4=RAk4.Jacobian*lambdaA4;
        J(nA+nB+10:nA+nB+11,1:nA)=tempa4(1:2,:);

        tempb4=RBk4.Jacobian*lambdaB4;
        J(nA+nB+10:nA+nB+11,nA+1:nA+nB)=-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*tempb4(1:2,:);
        J(nA+nB+10:nA+nB+11,nA+nB+11)=[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[-Lcon4*sin(gamma4);Lcon4*cos(gamma4)];


        % 受力平衡方程的导数
        % fa
        % fa对theta_a的导数
        % 求和的第0项 对应Fb
        tempA0=RA.partial;
        % 求和的第1项 对应fcon1
        tempA1=transpose(lambdaA1)*RAk1.partial*lambdaA1;
        % 求和的第2项 对应fcon2
        tempA2=transpose(lambdaA2)*RAk2.partial*lambdaA2;
        % 求和的第3项 对应fcon3
        tempA3=transpose(lambdaA3)*RAk3.partial*lambdaA3;
        % 求和的第4项 对应fcon4
        tempA4=transpose(lambdaA4)*RAk4.partial*lambdaA4;

        J(4:nA+3,1:nA)=RA.K_theta-1*(1*tempA0+tempA1+tempA2+tempA3+tempA4);

        % fa对F_b的导数
        J(4:nA+3,nA+nB+1:nA+nB+3)=-transpose(RA.Jacobian)*B;

        % fa对第1组约束的导数
        J(4:nA+3,nA+nB+4)=-transpose(lambdaA1)*transpose(RAk1.Jacobian)*(-rotz(-alpha_degree)*[cos(gamma1);sin(gamma1);0]);
        J(4:nA+3,nA+nB+5)=-transpose(lambdaA1)*transpose(RAk1.Jacobian)*(-rotz(-alpha_degree)*[-fcon1*sin(gamma1);fcon1*cos(gamma1);0]);
        % fa对第2组约束的导数
        J(4:nA+3,nA+nB+6)=-transpose(lambdaA2)*transpose(RAk2.Jacobian)*(-rotz(-alpha_degree)*[cos(gamma2);sin(gamma2);0]);
        J(4:nA+3,nA+nB+7)=-transpose(lambdaA2)*transpose(RAk2.Jacobian)*(-rotz(-alpha_degree)*[-fcon2*sin(gamma2);fcon2*cos(gamma2);0]);
        % fa对第3组约束的导数
        J(4:nA+3,nA+nB+8)=-transpose(lambdaA3)*transpose(RAk3.Jacobian)*(-rotz(-alpha_degree)*[cos(gamma3);sin(gamma3);0]);
        J(4:nA+3,nA+nB+9)=-transpose(lambdaA3)*transpose(RAk3.Jacobian)*(-rotz(-alpha_degree)*[-fcon3*sin(gamma3);fcon3*cos(gamma3);0]);
        % fa对第4组约束的导数
        J(4:nA+3,nA+nB+10)=-transpose(lambdaA4)*transpose(RAk4.Jacobian)*(-rotz(-alpha_degree)*[cos(gamma4);sin(gamma4);0]);
        J(4:nA+3,nA+nB+11)=-transpose(lambdaA4)*transpose(RAk4.Jacobian)*(-rotz(-alpha_degree)*[-fcon4*sin(gamma4);fcon4*cos(gamma4);0]);

        % fb
        % fb对theta_b的导数
        % 求和的第0项 对应Fb
        tempB0=RB.partial;
        % 求和的第1项 对应fcon1
        tempB1=transpose(lambdaB1)*RBk1.partial*lambdaB1;
        % 求和的第2项 对应fcon2
        tempB2=transpose(lambdaB2)*RBk2.partial*lambdaB2;
        % 求和的第3项 对应fcon3
        tempB3=transpose(lambdaB3)*RBk3.partial*lambdaB3;
        % 求和的第4项 对应fcon4
        tempB4=transpose(lambdaB4)*RBk4.partial*lambdaB4;

        J(nA+4:nA+nB+3,nA+1:nA+nB)=RB.K_theta-1*(1*tempB0+tempB1+tempB2+tempB3+tempB4);

        % fb对F_b的导数
        J(nA+4:nA+nB+3,nA+nB+1:nA+nB+3)=-transpose(RB.Jacobian);

        % fb对第1组约束的导数
        J(nA+4:nA+nB+3,nA+nB+4)=-transpose(lambdaB1)*transpose(RBk1.Jacobian)*(rotz(-beta_degree)*[cos(gamma1);sin(gamma1);0]);
        J(nA+4:nA+nB+3,nA+nB+5)=-transpose(lambdaB1)*transpose(RBk1.Jacobian)*(rotz(-beta_degree)*[-fcon1*sin(gamma1);fcon1*cos(gamma1);0]);
        % fb对第2组约束的导数
        J(nA+4:nA+nB+3,nA+nB+6)=-transpose(lambdaB2)*transpose(RBk2.Jacobian)*(rotz(-beta_degree)*[cos(gamma2);sin(gamma2);0]);
        J(nA+4:nA+nB+3,nA+nB+7)=-transpose(lambdaB2)*transpose(RBk2.Jacobian)*(rotz(-beta_degree)*[-fcon2*sin(gamma2);fcon2*cos(gamma2);0]);
        % fb对第3组约束的导数
        J(nA+4:nA+nB+3,nA+nB+8)=-transpose(lambdaB3)*transpose(RBk3.Jacobian)*(rotz(-beta_degree)*[cos(gamma3);sin(gamma3);0]);
        J(nA+4:nA+nB+3,nA+nB+9)=-transpose(lambdaB3)*transpose(RBk3.Jacobian)*(rotz(-beta_degree)*[-fcon3*sin(gamma3);fcon3*cos(gamma3);0]);
        % fb对第4组约束的导数
        J(nA+4:nA+nB+3,nA+nB+10)=-transpose(lambdaB4)*transpose(RBk4.Jacobian)*(rotz(-beta_degree)*[cos(gamma4);sin(gamma4);0]);
        J(nA+4:nA+nB+3,nA+nB+11)=-transpose(lambdaB4)*transpose(RBk4.Jacobian)*(rotz(-beta_degree)*[-fcon4*sin(gamma4);fcon4*cos(gamma4);0]);
    end
end