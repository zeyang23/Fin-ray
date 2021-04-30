% 真实finray的参数

% 单位约定
% 长度单位取mm
% 力的单位取kN
% 杨氏模量的单位取GPa

% 21-04-09
% 离谱的bug
% 右侧的刚性约束和柔性板的位置不协调

% 21-04-10 0:28
% 发现bug
% 右侧刚性约束的柔性板的长度设置错误 
% 原来是 kbi/obj.nB*obj.LA
% 现修改为 kbi*obj.RodB.seg_length

clear
clc


xA=0;
yA=0;

xB=73.24;
yB=0;


alpha_degree=90;
beta_degree=116;
psi_degree=beta_degree-alpha_degree;

psi=psi_degree/180*pi;
alpha=alpha_degree/180*pi;
beta=beta_degree/180*pi;


nA=50;
nB=50;


LA=(xB-xA)/tan(psi);
LB=(xB-xA)/sin(psi);

wid=28;
thi=0.2;
E=197;


constraint_ratio=[];
Lcon1=40.76;
Lcon2=35.78;
Lcon3=31.03;
Lcon4=26.37;
Lcon5=21.62;
Lcon6=16.6;

La_1=27;
La_2=54;
La_3=81;
La_4=108;

Lb_1=45.67;
Lb_2=78.69;
Lb_3=107.5;
Lb_4=132.27;

Lcon1=sqrt((LA-La_1)^2+(LB-Lb_1)^2-2*(LA-La_1)*(LB-Lb_1)*cos(psi));
Lcon2=sqrt((LA-La_2)^2+(LB-Lb_2)^2-2*(LA-La_2)*(LB-Lb_2)*cos(psi));
Lcon3=sqrt((LA-La_3)^2+(LB-Lb_3)^2-2*(LA-La_3)*(LB-Lb_3)*cos(psi));
Lcon4=sqrt((LA-La_4)^2+(LB-Lb_4)^2-2*(LA-La_4)*(LB-Lb_4)*cos(psi));

constraint_ratio=[Lcon1,La_1/LA,Lb_1/LB;
                  Lcon2,La_2/LA,Lb_2/LB;
                  Lcon3,La_3/LA,Lb_3/LB;
                  Lcon4,La_4/LA,Lb_4/LB];


A_force_ratio=[];
A_force_ratio=[0.00,1/2];


B_force_ratio=[];



% 根据上面设定的参数生成参数结构体
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