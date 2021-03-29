% ʹ��ga������������滮����

% 21-03-29
% ��1��������Գɹ��ˣ�exitFlagΪ1��˵��ga���ɹ�
% ga�ܹ�������������k��F
% ������������ʱ��ܳ�

% ��1������õ�������������Ϊ5�������F����Ϊ5����������Ϊ5�е�С
% ���������㷨ִ��ʱ��̫������ʱ�������µĳ��ԡ�

clear
clc

L0=1;

radius=0.1*L0;
center_x=0.2*L0;
center_y=0.5*L0;

f=@(tangent_var) myfunc(tangent_var,center_x,center_y,radius);

lb = [5,0];
ub = [45,5];

[X,fval,exitflag] = ga(f,2,[],[],[],[],lb,ub,[],1,[]);


% ��֤���
% ����Բ��λ��
rectangle('Position',[center_x-radius,center_y-radius,2*radius,2*radius],'Curvature',[1,1],'linewidth',1,'edgecolor','r')
axis equal;
hold on


tangent_ratio=X(1)/50;
tangent_F=X(2);

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


nA=50;
nB=50;


LA=1*L0;
LB=1*L0;

wid=5e-3;
thi=1e-3;
E=197*1e9;


constraint_ratio=[];
Lcon1=6/7*sqrt((xA-xB)^2+(yA-yB)^2);
Lcon2=5/7*sqrt((xA-xB)^2+(yA-yB)^2);
Lcon3=4/7*sqrt((xA-xB)^2+(yA-yB)^2);
Lcon4=3/7*sqrt((xA-xB)^2+(yA-yB)^2);
Lcon5=2/7*sqrt((xA-xB)^2+(yA-yB)^2);
Lcon6=1/7*sqrt((xA-xB)^2+(yA-yB)^2);
constraint_ratio=[Lcon1,1/7,1/7;
              Lcon2,2/7,2/7;
              Lcon3,3/7,3/7;
              Lcon4,4/7,4/7;
              Lcon5,5/7,5/7;
              Lcon6,6/7,6/7];


A_force_ratio=[];
A_force_ratio=[tangent_F,tangent_ratio];


B_force_ratio=[];


% ���������趨�Ĳ������ɲ����ṹ��
finray_info=struct();

finray_info.pA=[xA;yA;alpha];
finray_info.pB=[xB;yA;beta];

finray_info.LA=LA;
finray_info.LB=LB;

finray_info.nA=nA;
finray_info.nB=nB;

finray_info.psi=psi;

finray_info.wid_A=wid;
finray_info.wid_B=wid;

finray_info.thi_A=thi;
finray_info.thi_B=thi;

finray_info.E_A=E;
finray_info.E_B=E;

finray_info.constraint_ratio=constraint_ratio;

finray_info.A_force_ratio=A_force_ratio;

finray_info.B_force_ratio=B_force_ratio;


Finray1=finray(finray_info);




x0=zeros(nA+nB+3+2*Finray1.constraint_number,1);

options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

g=@(x) Finray1.cal_balance(x);
[x_solve,fval,exitflag,output] = fsolve(g,x0,options);

Finray1.plot_state(x_solve);




% ���������ݶȵ������ֱ�����Finray��Բ�����е�����

% ����Բ��λ��x y Բ���뾶r

% ���룺����λ��L ���Ա�������ʽ��������0��1֮�䣩 �е㴦�ķ�����F
% ������в�r1 r2
% ԭ��ʹ�в�Ϊ0��L��F������ʵ�����

function normr=myfunc(tangent_var,center_x,center_y,radius)

    tangent_index=tangent_var(1);
    
    tangent_F=tangent_var(2);

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


    nA=50;
    nB=50;


    LA=1*L0;
    LB=1*L0;

    wid=5e-3;
    thi=1e-3;
    E=197*1e9;


    constraint_ratio=[];
    Lcon1=6/7*sqrt((xA-xB)^2+(yA-yB)^2);
    Lcon2=5/7*sqrt((xA-xB)^2+(yA-yB)^2);
    Lcon3=4/7*sqrt((xA-xB)^2+(yA-yB)^2);
    Lcon4=3/7*sqrt((xA-xB)^2+(yA-yB)^2);
    Lcon5=2/7*sqrt((xA-xB)^2+(yA-yB)^2);
    Lcon6=1/7*sqrt((xA-xB)^2+(yA-yB)^2);
    constraint_ratio=[Lcon1,1/7,1/7;
                      Lcon2,2/7,2/7;
                      Lcon3,3/7,3/7;
                      Lcon4,4/7,4/7;
                      Lcon5,5/7,5/7;
                      Lcon6,6/7,6/7];


    A_force_ratio=[];
    tangent_ratio=tangent_index/nA;
    A_force_ratio=[tangent_F,tangent_ratio];


    B_force_ratio=[];


    % ���������趨�Ĳ������ɲ����ṹ��
    finray_info=struct();

    finray_info.pA=[xA;yA;alpha];
    finray_info.pB=[xB;yA;beta];

    finray_info.LA=LA;
    finray_info.LB=LB;

    finray_info.nA=nA;
    finray_info.nB=nB;

    finray_info.psi=psi;

    finray_info.wid_A=wid;
    finray_info.wid_B=wid;

    finray_info.thi_A=thi;
    finray_info.thi_B=thi;

    finray_info.E_A=E;
    finray_info.E_B=E;

    finray_info.constraint_ratio=constraint_ratio;

    finray_info.A_force_ratio=A_force_ratio;

    finray_info.B_force_ratio=B_force_ratio;


    Finray1=finray(finray_info);
    Finray1.A_force_number=size(Finray1.A_force_ratio,1);
    Finray1.A_force_index=zeros(size(Finray1.A_force_ratio));
    Finray1.A_force_index(:,1)=Finray1.A_force_ratio(:,1);


    Finray1.A_force_index(1,2)=tangent_index;
    Finray1.A_force_array(1)=planar_nR(Finray1.E_A,tangent_index/Finray1.nA*Finray1.LA,Finray1.wid_A,Finray1.thi_A,tangent_index,[0;0;0]);




    x0=zeros(nA+nB+3+2*Finray1.constraint_number,1);

    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

    f=@(x) Finray1.cal_balance(x);
    [x_solve,fval,exitflag,output] = fsolve(f,x0,options);
    
    
    
    % �õ��е��x y phi
    
    thetaA=x_solve(1:Finray1.nA);
    Finray1.A_force_array(1).theta=thetaA(1:Finray1.A_force_index(1,2));
    Finray1.A_force_array(1).cal_pe;

    pka=Finray1.A_force_array(1).pe;
    PA(1)=pka(1)*cos(Finray1.pA(3))-pka(2)*sin(Finray1.pA(3))+Finray1.pA(1);
    PA(2)=pka(1)*sin(Finray1.pA(3))+pka(2)*cos(Finray1.pA(3))+Finray1.pA(2);
    PA(3)=pka(3)+Finray1.pA(3);
    
    r=zeros(2,1);
    r(1)=center_x+radius*sin(PA(3))-PA(1);
    r(2)=center_y-radius*cos(PA(3))-PA(2);
    
    normr=norm(r);
end