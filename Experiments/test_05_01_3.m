% 5月1日实验
% Finray的形状重建

clear
clc

%% 标定第I块弯曲传感器
% 测试条件：电流4mA
volt1 = [];
volt1(1,:) = [100000000000000,2.08];
volt1(2,:) = [272,2.17];
volt1(3,:) = [203,2.21];
volt1(4,:) = [162,2.25];
volt1(5,:) = [134,2.305];
volt1(6,:) = [115,2.35];

volt1(7,:) = [100,2.40];
volt1(8,:) = [89,2.43];
volt1(9,:) = [80,2.49];
volt1(10,:) = [73,2.55];
volt1(11,:) = [67,2.58];



% 减掉初值
init1=volt1(1,2);
volt1(:,2)=volt1(:,2)-init1;


volt1(:,1) = 1./volt1(:,1)*38;
% plot(volt1(:,2),volt1(:,1),'or','linewidth',1);

coff_1 = polyfit(volt1(:,2),volt1(:,1),1);

% hold on
% grid on
% y_1 = polyval(coff_1,volt1(:,2));
% plot(volt1(:,2),y_1,'--r','linewidth',1)


%% 标定第II块弯曲传感器
volt2 = [];
volt2(1,:) = [100000000000000,1.81];
volt2(2,:) = [272,1.89];
volt2(3,:) = [203,1.92];
volt2(4,:) = [162,1.95];
volt2(5,:) = [134,1.98];
volt2(6,:) = [115,2.015];

volt2(7,:) = [100,2.05];
volt2(8,:) = [89,2.08];
volt2(9,:) = [80,2.12];
volt2(10,:) = [73,2.14];
volt2(11,:) = [67,2.17];


% 减掉初值
init2=volt2(1,2);
volt2(:,2)=volt2(:,2)-init2;


volt2(:,1) = 1./volt2(:,1)*38;
% plot(volt2(:,2),volt2(:,1),'og','linewidth',1);

coff_2 = polyfit(volt2(:,2),volt2(:,1),1);

% hold on
% grid on
% y_2 = polyval(coff_2,volt2(:,2));
% plot(volt2(:,2),y_2,'--g','linewidth',1)

%% 根据传感器示数，求出夹角

% % 初始测试
% % U1_0=2.01;
% % U2_0=1.56;
% % 第1组数据
% U1=2.07;
% U2=3.65-2.07;
% % 第2组数据
% U1=2.24;
% U2=1.69;


% 连续抓取测试
% U1_0=2.02;
% U2_0=1.58;
% % 第1组数据
% U1=2.07;
% U2=1.605;
% % 第2组数据
% U1=2.115;
% U2=1.63;
% % 第3组数据
% U1=2.14;
% U2=1.635;
% % 第4组数据
% U1=2.155;
% U2=1.645;
% % 第5组数据
% U1=2.19;
% U2=1.66;
% % 第6组数据
% U1=2.225;
% U2=1.69;
% % 第7组数据
% U1=2.247;
% U2=1.70;

% 回程
% % 第8组数据
% U1=2.21;
% U2=1.67;
% % 第9组数据
% U1=2.18;
% U2=1.66;
% % 第10组数据
% U1=2.14;
% U2=1.63;
% % 第11组数据
% U1=2.11;
% U2=1.62;
% % 第12组数据
% U1=2.04;
% U2=1.58;
% % 第13组数据
% U1=2.02;
% U2=1.58;


U1_delta=U1-U1_0;
U2_delta=U2-U2_0;


DeltaA=polyval(coff_1,U1_delta);
DeltaB=polyval(coff_2,U2_delta);


%% 反解Fin-ray形状
% 方法一：不需要触觉传感器信息

f=@(tangent_var) myfunc_1(tangent_var,DeltaA,DeltaB);

lb = [0,0];
ub = [1,8];
X0 = [0.5,1];
[X,resnorm,residual,exitflag_lsq,output_lsq] = lsqnonlin(f,X0,lb,ub);
% [X,fval_fsolve,exitflag_fsolve,output_fsolve]=fsolve(f,X0);



% 画出求解结果
tangent_ratio=X(1);
tangent_F=X(2);

xA=0;
yA=0;

xB=73.24e-3;
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

wid=28e-3;
thi=0.2e-3;
E=197e9;


constraint_ratio=[];

La_1=27e-3;
La_2=54e-3;
La_3=81e-3;
La_4=108e-3;

Lb_1=45.67e-3;
Lb_2=78.69e-3;
Lb_3=107.5e-3;
Lb_4=132.27e-3;

Lcon1=sqrt((LA-La_1)^2+(LB-Lb_1)^2-2*(LA-La_1)*(LB-Lb_1)*cos(psi));
Lcon2=sqrt((LA-La_2)^2+(LB-Lb_2)^2-2*(LA-La_2)*(LB-Lb_2)*cos(psi));
Lcon3=sqrt((LA-La_3)^2+(LB-Lb_3)^2-2*(LA-La_3)*(LB-Lb_3)*cos(psi));
Lcon4=sqrt((LA-La_4)^2+(LB-Lb_4)^2-2*(LA-La_4)*(LB-Lb_4)*cos(psi));

constraint_ratio=[Lcon1,La_1/LA,Lb_1/LB;
                  Lcon2,La_2/LA,Lb_2/LB;
                  Lcon3,La_3/LA,Lb_3/LB;
                  Lcon4,La_4/LA,Lb_4/LB];

A_force_ratio=[];
A_force_ratio=[tangent_F,tangent_ratio];


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


Finray_solve=finray_force(finray_info);




x0=zeros(nA+nB+3+2*Finray_solve.constraint_number,1);

options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

f=@(x) Finray_solve.cal_balance(x);
[x_solve,fval,exitflag,output] = fsolve(f,x0,options);

figure
Finray_solve.plot_state(x_solve);


% 传感器的长度
sensor_length=38e-3;
sensor_length_ratio=sensor_length/LB;

% 传感器的中心位置
sensor1_center_ratio=89.15e-3/LB;
sensor2_center_ratio=131.25e-3/LB;


p1=fix((sensor1_center_ratio-sensor_length_ratio/2)*nB);
q1=fix((sensor1_center_ratio+sensor_length_ratio/2)*nB);

p2=fix((sensor2_center_ratio-sensor_length_ratio/2)*nB);
q2=fix((sensor2_center_ratio+sensor_length_ratio/2)*nB);



Finray_solve.RodB.cal_posall;
sensor1_pos=Finray_solve.RodB.pos_all(1+p1:1+q1,:);
sensor2_pos=Finray_solve.RodB.pos_all(1+p2:1+q2,:);

plot_sensor_pos(sensor1_pos,beta,[xB,yB]);
plot_sensor_pos(sensor2_pos,beta,[xB,yB]);


%% 反解Fin-ray形状
% 方法二：需要触觉传感器信息
tangent_ratio_real=0;
f=@(tangent_var) myfunc_2(tangent_var,DeltaA,DeltaB,tangent_ratio_real);

lb = 0;
ub = 8;
X0 = 1;
[X,resnorm,residual,exitflag_lsq,output_lsq] = lsqnonlin(f,X0,lb,ub);
% [X,fval_fsolve,exitflag_fsolve,output_fsolve]=fsolve(f,X0);



% 求解结果
tangent_ratio=tangent_ratio_real;
tangent_F=X;

xA=0;
yA=0;

xB=73.24e-3;
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

wid=28e-3;
thi=0.2e-3;
E=197e9;


constraint_ratio=[];

La_1=27e-3;
La_2=54e-3;
La_3=81e-3;
La_4=108e-3;

Lb_1=45.67e-3;
Lb_2=78.69e-3;
Lb_3=107.5e-3;
Lb_4=132.27e-3;

Lcon1=sqrt((LA-La_1)^2+(LB-Lb_1)^2-2*(LA-La_1)*(LB-Lb_1)*cos(psi));
Lcon2=sqrt((LA-La_2)^2+(LB-Lb_2)^2-2*(LA-La_2)*(LB-Lb_2)*cos(psi));
Lcon3=sqrt((LA-La_3)^2+(LB-Lb_3)^2-2*(LA-La_3)*(LB-Lb_3)*cos(psi));
Lcon4=sqrt((LA-La_4)^2+(LB-Lb_4)^2-2*(LA-La_4)*(LB-Lb_4)*cos(psi));

constraint_ratio=[Lcon1,La_1/LA,Lb_1/LB;
                  Lcon2,La_2/LA,Lb_2/LB;
                  Lcon3,La_3/LA,Lb_3/LB;
                  Lcon4,La_4/LA,Lb_4/LB];

A_force_ratio=[];
A_force_ratio=[tangent_F,tangent_ratio];


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


Finray_solve=finray_force(finray_info);




x0=zeros(nA+nB+3+2*Finray_solve.constraint_number,1);

options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

f=@(x) Finray_solve.cal_balance(x);
[x_solve,fval,exitflag,output] = fsolve(f,x0,options);

figure
Finray_solve.plot_state(x_solve);


% 传感器的长度
sensor_length=38e-3;
sensor_length_ratio=sensor_length/LB;

% 传感器的中心位置
sensor1_center_ratio=90e-3/LB;
sensor2_center_ratio=135e-3/LB;


p1=fix((sensor1_center_ratio-sensor_length_ratio/2)*nB);
q1=fix((sensor1_center_ratio+sensor_length_ratio/2)*nB);

p2=fix((sensor2_center_ratio-sensor_length_ratio/2)*nB);
q2=fix((sensor2_center_ratio+sensor_length_ratio/2)*nB);



Finray_solve.RodB.cal_posall;
sensor1_pos=Finray_solve.RodB.pos_all(1+p1:1+q1,:);
sensor2_pos=Finray_solve.RodB.pos_all(1+p2:1+q2,:);

plot_sensor_pos(sensor1_pos,beta,[xB,yB]);
plot_sensor_pos(sensor2_pos,beta,[xB,yB]);


%%
function r=myfunc_2(tangent_var,Delta1,Delta2,tangent_ratio)

    tangent_F=tangent_var(1);


    xA=0;
    yA=0;

    xB=73.24e-3;
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

    wid=28e-3;
    thi=0.2e-3;
    E=197e9;


    constraint_ratio=[];

    La_1=27e-3;
    La_2=54e-3;
    La_3=81e-3;
    La_4=108e-3;

    Lb_1=45.67e-3;
    Lb_2=78.69e-3;
    Lb_3=107.5e-3;
    Lb_4=132.27e-3;

    Lcon1=sqrt((LA-La_1)^2+(LB-Lb_1)^2-2*(LA-La_1)*(LB-Lb_1)*cos(psi));
    Lcon2=sqrt((LA-La_2)^2+(LB-Lb_2)^2-2*(LA-La_2)*(LB-Lb_2)*cos(psi));
    Lcon3=sqrt((LA-La_3)^2+(LB-Lb_3)^2-2*(LA-La_3)*(LB-Lb_3)*cos(psi));
    Lcon4=sqrt((LA-La_4)^2+(LB-Lb_4)^2-2*(LA-La_4)*(LB-Lb_4)*cos(psi));

    constraint_ratio=[Lcon1,La_1/LA,Lb_1/LB;
                      Lcon2,La_2/LA,Lb_2/LB;
                      Lcon3,La_3/LA,Lb_3/LB;
                      Lcon4,La_4/LA,Lb_4/LB];

    A_force_ratio=[];
    A_force_ratio=[tangent_F,tangent_ratio];


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


    FinrayA=finray_force(finray_info);




    x0=zeros(nA+nB+3+2*FinrayA.constraint_number,1);

    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

    f=@(x) FinrayA.cal_balance(x);
    [x_solve,fval,exitflag,output] = fsolve(f,x0,options);
    
    
    % 传感器的长度
    sensor_length=38e-3;
    sensor_length_ratio=sensor_length/LB;

    % 传感器的中心位置
    sensor1_center_ratio=90e-3/LB;
    sensor2_center_ratio=135e-3/LB;


    p1=fix((sensor1_center_ratio-sensor_length_ratio/2)*nB);
    q1=fix((sensor1_center_ratio+sensor_length_ratio/2)*nB);

    p2=fix((sensor2_center_ratio-sensor_length_ratio/2)*nB);
    q2=fix((sensor2_center_ratio+sensor_length_ratio/2)*nB);


    Delta1_now=sum(FinrayA.RodB.theta(p1:q1));
    Delta2_now=sum(FinrayA.RodB.theta(p2:q2));

    
    
    
    % 拿到切点的x y phi
    
    
    r=zeros(2,1);
    r(1)=Delta1_now-Delta1;
    r(2)=Delta2_now-Delta2;
    
end



%%
function r=myfunc_1(tangent_var,Delta1,Delta2)

    tangent_ratio=tangent_var(1);
    tangent_F=tangent_var(2);


    xA=0;
    yA=0;

    xB=73.24e-3;
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

    wid=28e-3;
    thi=0.2e-3;
    E=197e9;


    constraint_ratio=[];

    La_1=27e-3;
    La_2=54e-3;
    La_3=81e-3;
    La_4=108e-3;

    Lb_1=45.67e-3;
    Lb_2=78.69e-3;
    Lb_3=107.5e-3;
    Lb_4=132.27e-3;

    Lcon1=sqrt((LA-La_1)^2+(LB-Lb_1)^2-2*(LA-La_1)*(LB-Lb_1)*cos(psi));
    Lcon2=sqrt((LA-La_2)^2+(LB-Lb_2)^2-2*(LA-La_2)*(LB-Lb_2)*cos(psi));
    Lcon3=sqrt((LA-La_3)^2+(LB-Lb_3)^2-2*(LA-La_3)*(LB-Lb_3)*cos(psi));
    Lcon4=sqrt((LA-La_4)^2+(LB-Lb_4)^2-2*(LA-La_4)*(LB-Lb_4)*cos(psi));

    constraint_ratio=[Lcon1,La_1/LA,Lb_1/LB;
                      Lcon2,La_2/LA,Lb_2/LB;
                      Lcon3,La_3/LA,Lb_3/LB;
                      Lcon4,La_4/LA,Lb_4/LB];

    A_force_ratio=[];
    A_force_ratio=[tangent_F,tangent_ratio];


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


    FinrayA=finray_force(finray_info);




    x0=zeros(nA+nB+3+2*FinrayA.constraint_number,1);

    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

    f=@(x) FinrayA.cal_balance(x);
    [x_solve,fval,exitflag,output] = fsolve(f,x0,options);
    
    
    % 传感器的长度
    sensor_length=38e-3;
    sensor_length_ratio=sensor_length/LB;

    % 传感器的中心位置
    sensor1_center_ratio=90e-3/LB;
    sensor2_center_ratio=135e-3/LB;


    p1=fix((sensor1_center_ratio-sensor_length_ratio/2)*nB);
    q1=fix((sensor1_center_ratio+sensor_length_ratio/2)*nB);

    p2=fix((sensor2_center_ratio-sensor_length_ratio/2)*nB);
    q2=fix((sensor2_center_ratio+sensor_length_ratio/2)*nB);


    Delta1_now=sum(FinrayA.RodB.theta(p1:q1));
    Delta2_now=sum(FinrayA.RodB.theta(p2:q2));

    
    
    
    % 拿到切点的x y phi
    
    
    r=zeros(2,1);
    r(1)=Delta1_now-Delta1;
    r(2)=Delta2_now-Delta2;
    
end


function abs_pos=plot_sensor_pos(pos_all,alpha,pA)
    xA=pA(1);
    yA=pA(2);
    abs_pos=zeros(size(pos_all));
    n=size(pos_all,1);
    for i=1:n
        abs_pos(i,1)=pos_all(i,1)*cos(alpha)-pos_all(i,2)*sin(alpha)+xA;
        abs_pos(i,2)=pos_all(i,1)*sin(alpha)+pos_all(i,2)*cos(alpha)+yA;
    end
    plot(abs_pos(:,1),abs_pos(:,2),'-o','MarkerSize',6);
    axis equal
end