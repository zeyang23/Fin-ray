% ��ʵfinray�Ĳ���

% ��λԼ��
% ���ȵ�λȡmm
% ���ĵ�λȡkN
% ����ģ���ĵ�λȡGPa

% 21-04-09
% ���׵�bug
% �Ҳ�ĸ���Լ�������԰��λ�ò�Э��

% 21-04-10 0:28
% ����bug
% �Ҳ����Լ�������԰�ĳ������ô��� 
% ԭ���� kbi/obj.nB*obj.LA
% ���޸�Ϊ kbi*obj.RodB.seg_length

clear
clc


xA=0;
yA=0;

xB=43.34;
yB=16.18;


psi_degree=21.29;
alpha_degree=90;
beta_degree=110.88;

psi=psi_degree/180*pi;
alpha=alpha_degree/180*pi;
beta=beta_degree/180*pi;


nA=50;
nB=50;


LA=129.32;
LB=121.16;

wid=43.4;
thi=5;

E=0.0094;


constraint_ratio=[];
Lcon1=40.76;
Lcon2=35.78;
Lcon3=31.03;
Lcon4=26.37;
Lcon5=21.62;
Lcon6=16.6;
constraint_ratio=[Lcon1,15.54/LA,16.71/LB;
                  Lcon2,30.54/LA,34.22/LB;
                  Lcon3,45.54/LA,50.71/LB;
                  Lcon4,60.54/LA,66.18/LB;
                  Lcon5,75.54/LA,80.6/LB;
                  Lcon6,90.54/LA,93.87/LB];


A_force_ratio=[];
A_force_ratio=[0,1/2;
               0,1/4];

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

options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false);

f=@(x) Finray1.cal_balance(x);
[x_solve,fval,exitflag,output] = fsolve(f,x0,options);

Finray1.plot_state(x_solve);