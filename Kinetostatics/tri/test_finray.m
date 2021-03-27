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


nA=50;
nB=50;


LA=1*L0;
LB=1*L0;

wid=5e-3;
thi=1e-3;
E=197*1e9;


Lcon1=4/5*sqrt((xA-xB)^2+(yA-yB)^2);
constraint_ratio=[Lcon1,1/5,1/5];
% constraint_ratio=[];

A_force_ratio=[1.5,1/2];
B_force_ratio=[2.5,1/4;0.5,3/4];



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

x0=zeros(nA+nB+5,1);

options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',true);

f=@(x) Finray1.cal_balance(x);
[x_solve,fval,exitflag,output] = fsolve(f,x0,options);

Finray1.plot_state(x_solve);

% Finray1.RodA.cal_posall;
% Finray1.RodB.cal_posall;
% 
% A_pos_all=Finray1.RodA.pos_all;
% B_pos_all=Finray1.RodB.pos_all;
% 
% A_abs_pos_all=plot_abs_pos(A_pos_all,alpha,[xA,yA]);
% hold on
% B_abs_pos_all=plot_abs_pos(B_pos_all,beta,[xB,yB]);
% 
% pka1=Finray1.A_constraint_array(1).pe;
% PA1(1)=pka1(1)*cos(alpha)-pka1(2)*sin(alpha)+xA;
% PA1(2)=pka1(1)*sin(alpha)+pka1(2)*cos(alpha)+yA;
% 
% pkb1=Finray1.B_constraint_array(1).pe;
% PB1(1)=pkb1(1)*cos(beta)-pkb1(2)*sin(beta)+xB;
% PB1(2)=pkb1(1)*sin(beta)+pkb1(2)*cos(beta)+yB;
% 
% plot([PA1(1) PB1(1)],[PA1(2) PB1(2)])