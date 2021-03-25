% ʹ��fsolve���
% �����ڲ���2������Լ��

% 21-03-24
% GradientCheck��Ȼ��ͨ��������
% ����fsolve�ṩ�ݶ�Ҳ�ܽ������ֻ�����ٶȺ���

% 21-03-24��
% �ƣ�֮ǰд����B=-rotz(beta-alpha)�����˰ѻ�����ת���ɽǶ����ˡ�����

% 21-03-25
% ����֪��GradientCheckΪɶͨ�����ˣ���Ȼ�����۹�ʽ�ƴ��ˡ�����
% ���ۣ�LZY̫����
% ���깫ʽ�Ժ�������ͨ��GradientCheck�ˡ�����
% �ܽ᣺������Ҫ���ͨ��Ҫ�����ӣ��������ӵĺ�������˷Ѵ���ʱ��

clear
clc

%% ��������
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


% ������ѡȡ
LA=0.8*L0;

% ��1�������Ĳ���
% ��1������Լ��
ka1=fix(1/8*nA);
kb1=fix(1/8*nB);
Lcon1=7/8*sqrt((xA-xB)^2+(yA-yB)^2);

LambdaA=zeros(nA,ka1);
LambdaA(1:ka1,1:ka1)=diag(ones(ka1,1));
LambdaB=zeros(nB,kb1);
LambdaB(1:kb1,1:kb1)=diag(ones(kb1,1));

wid=5e-3;
thi=1e-3;
E=197*1e9;

pdes=[0;0;0];

LB=L0;

RA=planar_nR(E,LA,wid,thi,nA,pdes);
RAk1=planar_nR(E,LA*ka1/nA,wid,thi,ka1,pdes);
RB=planar_nR(E,LB,wid,thi,nB,pdes);
RBk1=planar_nR(E,LB*kb1/nB,wid,thi,kb1,pdes);



%% fsolve��ⷽ����

x0=zeros(nA+nB+5,1);

%����ʹ���޸���Լ��ʱ�ĳ�ֵ
% load('x_init.mat')
% x0(1:nA+nB+3)=x_init;

options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',true);
% options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',true,'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-10);

f=@(x) myfun(x,xA,xB,yA,yB,beta,beta_degree,alpha,alpha_degree,A,B,b,nA,nB,RA,RB,RAk1,RBk1,Lcon1,ka1,kb1,LambdaA,LambdaB);
[x,fval,exitflag,output] = fsolve(f,x0,options);


%% �����ڲ��ĸ���Լ��
RA.cal_posall;
RB.cal_posall;

A_pos_all=RA.pos_all;
B_pos_all=RB.pos_all;

A_abs_pos_all=plot_abs_pos(A_pos_all,alpha,[xA,yA]);
hold on
B_abs_pos_all=plot_abs_pos(B_pos_all,beta,[xB,yB]);


% ��1������Լ��
pka1=RAk1.pe;
PA1(1)=pka1(1)*cos(alpha)-pka1(2)*sin(alpha)+xA;
PA1(2)=pka1(1)*sin(alpha)+pka1(2)*cos(alpha)+yA;

pkb1=RBk1.pe;
PB1(1)=pkb1(1)*cos(beta)-pkb1(2)*sin(beta)+xB;
PB1(2)=pkb1(1)*sin(beta)+pkb1(2)*cos(beta)+yB;

plot([PA1(1) PB1(1)],[PA1(2) PB1(2)])


%%
function [r,J]=myfun(x,xA,xB,yA,yB,beta,beta_degree,alpha,alpha_degree,A,B,b,nA,nB,RA,RB,RAk1,RBk1,Lcon,k1,k2,LambdaA,LambdaB)
    
    thetaA=zeros(nA,1);
    thetaKA=zeros(k1,1);
    thetaB=zeros(nB,1);
    thetaKB=zeros(k2,1);
    FB=zeros(3,1);
    fcon=0;
    gamma=0;
    
    thetaA=x(1:nA);
    thetaB=x(nA+1:nA+nB);
    FB=x(nA+nB+1:nA+nB+3);
    fcon=x(end-1);
    gamma=x(end);
    
    
    thetaKA=thetaA(1:k1);
    thetaKB=thetaB(1:k2);
    
    % ��ʼ���㺯��ֵ
    
    FA=B*FB;
    
    RA.theta=thetaA;
    RA.F=FA;
    RB.theta=thetaB;
    RB.F=FB;
    
    RA.update;
    RB.update;
    
    
    RAk1.theta=thetaKA;
    RBk1.theta=thetaKB;
    RAk1.F=(-rotz(-alpha_degree))*[fcon*cos(gamma);fcon*sin(gamma);0];
    RBk1.F=rotz(-beta_degree)*[fcon*cos(gamma);fcon*sin(gamma);0];
    RAk1.update;
    RBk1.update;
    
    pk1=RAk1.pe;
    pk2=RBk1.pe;
    
    
    r2=[pk1(1);pk1(2)]-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*[pk2(1);pk2(2)]+...
        [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[xA-xB+Lcon*cos(gamma);yA-yB+Lcon*sin(gamma)];
    
    r1=RA.pe-A*RB.pe-b;
    r3=RA.K_theta*thetaA-transpose(RA.Jacobian)*FA-LambdaA*transpose(RAk1.Jacobian)*RAk1.F;
    r4=RB.K_theta*thetaB-transpose(RB.Jacobian)*FB-LambdaB*transpose(RBk1.Jacobian)*RBk1.F;
    r=[r1;r2;r3;r4];
    
    
    % ��ʼ�����ſɱȾ���
    if nargout > 1
        J=zeros(nA+nB+5);

        J(1:3,1:nA)=RA.Jacobian;
        J(1:3,nA+1:nA+nB)=-A*RB.Jacobian;
        temp1=RAk1.Jacobian;
        J(4:5,1:k1)=temp1(1:2,:);
        
        temp2=RBk1.Jacobian;
        J(4:5,nA+1:nA+k2)=-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*temp2(1:2,:);
        J(4:5,end)=[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[-Lcon*sin(gamma);Lcon*cos(gamma)];

        temp3=RA.partial;
        temp4=LambdaA*RAk1.partial;
        J(6:nA+5,1:nA)=RA.K_theta-1*(1*temp3+temp4*transpose(LambdaA));

        J(6:nA+5,nA+nB+1:nA+nB+3)=-transpose(RA.Jacobian)*B;
        J(6:nA+5,end-1)=-LambdaA*transpose(RAk1.Jacobian)*(-rotz(-alpha_degree)*[cos(gamma);sin(gamma);0]);
        J(6:nA+5,end)=-LambdaB*transpose(RAk1.Jacobian)*(-rotz(-alpha_degree)*[-fcon*sin(gamma);fcon*cos(gamma);0]);

        temp5=RB.partial;
        temp6=LambdaB*RBk1.partial;
        J(nA+6:end,nA+1:nA+nB)=RB.K_theta-1*(1*temp5+temp6*transpose(LambdaB));

        J(nA+6:end,nA+nB+1:nA+nB+3)=-transpose(RB.Jacobian);
        J(nA+6:end,end-1)=-LambdaA*transpose(RBk1.Jacobian)*(rotz(-beta_degree)*[cos(gamma);sin(gamma);0]);
        J(nA+6:end,end)=-LambdaB*transpose(RBk1.Jacobian)*(rotz(-beta_degree)*[-fcon*sin(gamma);fcon*cos(gamma);0]);
    end
end