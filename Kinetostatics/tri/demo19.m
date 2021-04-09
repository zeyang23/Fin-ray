% 尝试建立kc与rc的单值函数关系

% 21-03-30
% 失败
% kc与rc一一对应的这个想法不太合理
% 如果目标切点的位置太靠下，这是不可能解的出来的

% 梯度检测也没有通过，不清楚为什么

% 21-03-30中午
% 放弃这个思路


clear
clc

L0=1;

xc=0.2*L0;
yc=0.5*L0;

rc=0.1*L0;

A_contact_index=10;


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


nA=20;
nB=20;


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
% constraint_ratio=[Lcon1,1/7,1/7;
%                   Lcon2,2/7,2/7;
%                   Lcon3,3/7,3/7;
%                   Lcon4,4/7,4/7;
%                   Lcon5,5/7,5/7;
%                   Lcon6,6/7,6/7];

constraint_ratio=[Lcon3,3/7,3/7];



% 根据上面设定的参数生成参数结构体
finray_info=struct();

finray_info.xc=xc;
finray_info.yc=yc;

finray_info.rc=rc;

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

finray_info.A_contact_index=A_contact_index;

finray_info.A_contact_force=0;


Finray1=finray_contact(finray_info);


x0=zeros(nA+nB+5+2*Finray1.constraint_number,1);
x0(end)=rc;

options = optimoptions('fsolve','SpecifyObjectiveGradient',false,'CheckGradient',true);

f=@(x) Finray1.cal_balance(x);
[x_solve,fval,exitflag,output] = fsolve(f,x0,options);

Finray1.plot_state(x_solve);

rectangle('Position',[xc-rc,yc-rc,2*rc,2*rc],'Curvature',[1,1],'linewidth',1,'edgecolor','r')
axis equal;