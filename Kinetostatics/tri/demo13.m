% ʹ��fsolve���
% �����ڲ���2������Լ��

% 21-03-24
% GradientCheck��Ȼ��ͨ��������
% ����fsolve�ṩ�ݶ�Ҳ�ܽ������ֻ�����ٶȺ���

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

B=-rotz(beta-alpha);


nA=50;
nB=50;


% ������ѡȡ
LA=0.8*L0;

% ��1�������Ĳ���
% ��1������Լ��
ka1=fix(1/8*nA);
kb1=fix(1/8*nB);
Lcon1=7/8*sqrt((xA-xB)^2+(yA-yB)^2);

vA1=[ones(ka1,1);zeros(nA-ka1,1)];
lambdaA1=diag(vA1);
vB1=[ones(kb1,1);zeros(nB-kb1,1)];
lambdaB1=diag(vB1);
       


wid=5e-3;
thi=1e-3;
E=197*1e9;

pdes=[0;0;0];

LB=L0;

RA=planar_nR(E,LA,wid,thi,nA,pdes);
RB=planar_nR(E,LB,wid,thi,nB,pdes);




%% fsolve��ⷽ����

x0=zeros(nA+nB+5,1);

%����ʹ���޸���Լ��ʱ�ĳ�ֵ
load('x_init.mat')
x0(1:nA+nB+3)=x_init;

options = optimoptions('fsolve','SpecifyObjectiveGradient',false,'CheckGradient',false);

f=@(x) myfun(x,xA,xB,yA,yB,beta,beta_degree,alpha,alpha_degree,A,B,b,nA,nB,RA,RB,Lcon1,ka1,kb1,lambdaA1,lambdaB1);
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
pka1=RA.cal_pk(ka1);
PA1(1)=pka1(1)*cos(alpha)-pka1(2)*sin(alpha)+xA;
PA1(2)=pka1(1)*sin(alpha)+pka1(2)*cos(alpha)+yA;

pkb1=RB.cal_pk(kb1);
PB1(1)=pkb1(1)*cos(beta)-pkb1(2)*sin(beta)+xB;
PB1(2)=pkb1(1)*sin(beta)+pkb1(2)*cos(beta)+yB;

plot([PA1(1) PB1(1)],[PA1(2) PB1(2)])


%%
function [r,J]=myfun(x,xA,xB,yA,yB,beta,beta_degree,alpha,alpha_degree,A,B,b,nA,nB,RA,RB,Lcon,k1,k2,lambdaA,lambdaB)
    
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
    
    
    % ��ʼ���㺯��ֵ
    
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
    
    
    % ��ʼ�����ſɱȾ���
    if nargout > 1
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
    end
end