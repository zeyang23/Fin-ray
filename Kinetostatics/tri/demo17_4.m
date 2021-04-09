% ����lsqnonlin����������������

clear
clc

L0=1;

radius=0.1*L0;
center_x=0.2*L0;
center_y=0.5*L0;

lb = [0,0];
ub = [1,10];

X0 = [0.5,1];

delta=linspace(0.1*L0,0,21);

X_series_lsq=zeros(21,3);

for i=1:21

    f=@(tangent_var) myfunc(tangent_var,delta(i),center_x,center_y,radius);

    

    [X,resnorm,residual,exitflag_lsq,output_lsq] = lsqnonlin(f,X0,lb,ub);
    
    X_series_lsq(i,:)=[X(1) X(2) resnorm];
    
    X0=X;
end

save('X_series_lsq.mat','X_series_lsq')




% ���������ݶȵ������ֱ�����Finray��Բ�����е�����

% ����Բ��λ��x y Բ���뾶r

% ���룺����λ��L ���Ա�������ʽ��������0��1֮�䣩 �е㴦�ķ�����F
% ������в�r1 r2
% ԭ��ʹ�в�Ϊ0��L��F������ʵ�����

function r=myfunc(tangent_var,delta,center_x,center_y,radius)

    tangent_ratio=tangent_var(1);
    tangent_F=tangent_var(2);

    L0=1;

    xA=0+delta;
    yA=0;

    xB=0.35*L0+delta;
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
    finray_info.pB=[xB;yB;beta];

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


    Finray1=finray_force(finray_info);




    x0=zeros(nA+nB+3+2*Finray1.constraint_number,1);

    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

    f=@(x) Finray1.cal_balance(x);
    [x_solve,fval,exitflag,output] = fsolve(f,x0,options);
    
    
    
    % �õ��е��x y phi
    
    thetaA=x_solve(1:Finray1.nA);
    Finray1.A_force_array(1).theta=thetaA(1:Finray1.A_force_index(1,2));
    Finray1.A_force_array(1).cal_pe;

    pka=Finray1.A_force_array(1).pe;
    L_tail=Finray1.A_force_ratio(1,2)*Finray1.LA-(Finray1.A_force_array(1).Ltotal-Finray1.A_force_array(1).seg_length/2);
    pka(1)=pka(1)-Finray1.A_force_array(1).seg_length*cos(sum(Finray1.A_force_array(1).theta))/2+L_tail*cos(sum(Finray1.A_force_array(1).theta));
    pka(2)=pka(2)-Finray1.A_force_array(1).seg_length*sin(sum(Finray1.A_force_array(1).theta))/2+L_tail*sin(sum(Finray1.A_force_array(1).theta));
    
    PA(1)=pka(1)*cos(Finray1.pA(3))-pka(2)*sin(Finray1.pA(3))+Finray1.pA(1);
    PA(2)=pka(1)*sin(Finray1.pA(3))+pka(2)*cos(Finray1.pA(3))+Finray1.pA(2);
    PA(3)=pka(3)+Finray1.pA(3);
    
    r=zeros(2,1);
    r(1)=center_x+radius*sin(PA(3))-PA(1);
    r(2)=center_y-radius*cos(PA(3))-PA(2);
    
end