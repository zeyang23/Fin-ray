% 测试多点接触 多点均相切

% 设置均匀间隔，只搜索起点和终点


%%
clear
clc

L0=1;

radius=0.2*L0;
center_x=0.1*L0;
center_y=0.4*L0;

f=@(tangent_var) myfunc(tangent_var,center_x,center_y,radius);


% 相切点的个数
N_tan=5;


lb=zeros(2+N_tan,1);

ub=zeros(N_tan+2,1);
ub(1)=1;
ub(2)=1;
for i=1:N_tan
    ub(i+2)=10;
end

ratio_series=linspace(0.2,0.8,N_tan);

X0=zeros(2+N_tan,1);
X0(1)=0.3;
X0(2)=0.7;
for i=1:N_tan
    X0(i+2)=1;
end


[X,resnorm,residual,exitflag_lsq,output_lsq] = lsqnonlin(f,X0,lb,ub);
% [X,fval_fsolve,exitflag_fsolve,output_fsolve]=fsolve(f,X0);


%% 验证结果
% 画出圆柱位置
rectangle('Position',[center_x-radius,center_y-radius,2*radius,2*radius],'Curvature',[1,1],'linewidth',1,'edgecolor','r')
axis equal;
hold on


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

start_ratio=X(1);
end_ratio=X(2);

N_tan=length(X)-2;

tangent_ratio=linspace(start_ratio,end_ratio,N_tan);
tangent_ratio=transpose(tangent_ratio);

tangent_F=zeros(N_tan,1);

for i=1:N_tan
    tangent_F(i)=X(i+2);
end

A_force_ratio=[tangent_F,tangent_ratio];


B_force_ratio=[];


% 根据上面设定的参数生成参数结构体
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


Finray1=finray_force(finray_info);




x0=zeros(nA+nB+3+2*Finray1.constraint_number,1);

options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

f=@(x) Finray1.cal_balance(x);
[x_solve,fval,exitflag,output] = fsolve(f,x0,options);

Finray1.plot_state(x_solve);




% 尝试在无梯度的情况下直接求解Finray与圆柱相切的问题

% 给定圆心位置x y 圆柱半径r

% 输入：相切位置L （以比例的形式给出，在0和1之间） 切点处的法向力F
% 输出：残差r1 r2
% 原理：使残差为0的L和F就是真实的情况

function r=myfunc(tangent_var,center_x,center_y,radius)

    start_ratio=tangent_var(1);
    end_ratio=tangent_var(2);
    
    N_tan=length(tangent_var)-2;
    
    tangent_ratio=linspace(start_ratio,end_ratio,N_tan);
    tangent_ratio=transpose(tangent_ratio);
    
    tangent_F=zeros(N_tan,1);
    
    for i=1:N_tan
        tangent_F(i)=tangent_var(i+2);
    end


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


    % 根据上面设定的参数生成参数结构体
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


    Finray1=finray_force(finray_info);




    x0=zeros(nA+nB+3+2*Finray1.constraint_number,1);

    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

    f=@(x) Finray1.cal_balance(x);
    [x_solve,fval,exitflag,output] = fsolve(f,x0,options);
    
    
    
    % 拿到切点的x y phi
    % 验证相切条件
    
    thetaA=x_solve(1:Finray1.nA);
    
    r=zeros(2*N_tan,1);
    
    for i=1:N_tan
        Finray1.A_force_array(i).theta=thetaA(1:Finray1.A_force_index(i,2));
        Finray1.A_force_array(i).cal_pe;

        pka=Finray1.A_force_array(i).pe;
        L_tail=Finray1.A_force_ratio(i,2)*Finray1.LA-(Finray1.A_force_array(i).Ltotal-Finray1.A_force_array(i).seg_length/2);
        pka(1)=pka(1)-Finray1.A_force_array(i).seg_length*cos(sum(Finray1.A_force_array(i).theta))/2+L_tail*cos(sum(Finray1.A_force_array(i).theta));
        pka(2)=pka(2)-Finray1.A_force_array(i).seg_length*sin(sum(Finray1.A_force_array(i).theta))/2+L_tail*sin(sum(Finray1.A_force_array(i).theta));

        PA(1)=pka(1)*cos(Finray1.pA(3))-pka(2)*sin(Finray1.pA(3))+Finray1.pA(1);
        PA(2)=pka(1)*sin(Finray1.pA(3))+pka(2)*cos(Finray1.pA(3))+Finray1.pA(2);
        PA(3)=pka(3)+Finray1.pA(3);
        
        r(2*i-1)=center_x+radius*sin(PA(3))-PA(1);
        r(2*i)=center_y-radius*cos(PA(3))-PA(2);
    end
    
end