clear
clc

% 形状重建
% 不规则物体

% 初始位置在旋转中心系中的坐标
% XC=30e-3;
% YC=-62.35e-3;

%% 根据传感器示数，求出夹角

coff_1=[1.186256198447459,0.056997358385406]; % 实际计算时取1.5
coff_2=[1.151462241145653,0.027617725599564]; % 实际计算时取1.5


U1_0=3.17-1.47;
U2_0=1.47;

% 第1组数据 0°delta=20
% U1=3.459-1.636;
% U2=1.636;

% 第2组数据 20°delta=20
% U1=3.47-1.64;
% U2=1.64;

% 第3组数据 40°delta=20
% U1=3.461-1.624;
% U2=1.624;

% 第4组数据 60°delta=25
% U1=3.536-1.66;
% U2=1.66;

% 第5组数据 80°delta=25
% U1=3.49-1.652;
% U2=1.652;

% 第6组数据 100°delta=25
% U1=3.51-1.676;
% U2=1.676;

% 第7组数据 120°delta=20
% U1=3.44-1.632;
% U2=1.632;

% 第8组数据 140°delta=20
% U1=3.461-1.64;
% U2=1.64;

% 第9组数据 160°delta=20
% U1=3.46-1.63;
% U2=1.63;

% 第10组数据 180°delta=20
% U1=3.433-1.62;
% U2=1.62;

% 第11组数据 200°delta=20
% U1=3.423-1.624;
% U2=1.624;

% 第12组数据 220°delta=20
% U1=3.45-1.64;
% U2=1.64;

% 第13组数据 240°delta=20
% U1=3.48-1.648;
% U2=1.648;

% 第14组数据 260°delta=20
% U1=3.48-1.644;
% U2=1.644;

% 第15组数据 280°delta=20
% U1=3.455-1.62;
% U2=1.62;

% 第16组数据 300°delta=25
% U1=3.51-1.649;
% U2=1.649;

% 第17组数据 320°delta=25
% U1=3.476-1.656;
% U2=1.656;

% 第18组数据 340°delta=20
% U1=3.42-1.627;
% U2=1.627;

% 第19组数据 360°delta=20
% U1=3.45-1.64;
% U2=1.64;


U1_delta=U1-U1_0;
U2_delta=U2-U2_0;


DeltaA=polyval(coff_1,U1_delta);
DeltaB=polyval(coff_2,U2_delta);


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
E=213e9;


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

finray_info.E_A=LA/135e-3*E;
finray_info.E_B=LB/156e-3*E;

finray_info.constraint_ratio=constraint_ratio;

finray_info.A_force_ratio=A_force_ratio;

finray_info.B_force_ratio=B_force_ratio;


Finray_solve=finray_force(finray_info);




x0=zeros(nA+nB+3+2*Finray_solve.constraint_number,1);

options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

f=@(x) Finray_solve.cal_balance(x);
[x_solve,fval,exitflag,output] = fsolve(f,x0,options);

figure
Finray_solve.plot_state_sr(x_solve);


% 传感器的长度
sensor_length=36e-3;

% 传感器的中心位置
sensor1_center=133.9e-3;
sensor2_center=93.1e-3;


p1=fix((sensor1_center-sensor_length/2-LB/nB/2)/(LB/nB))+1;
q1=fix((sensor1_center+sensor_length/2-LB/nB/2)/(LB/nB))+1;

p2=fix((sensor2_center-sensor_length/2-LB/nB/2)/(LB/nB))+1;
q2=fix((sensor2_center+sensor_length/2-LB/nB/2)/(LB/nB))+1;



Finray_solve.RodB.cal_posall;
sensor1_pos=Finray_solve.RodB.pos_all(1+p1:1+q1,:);
sensor2_pos=Finray_solve.RodB.pos_all(1+p2:1+q2,:);

plot_sensor_pos(sensor1_pos,beta,[xB,yB]);
plot_sensor_pos(sensor2_pos,beta,[xB,yB]);


axis equal
axis([-80e-3 100e-3 0 160e-3]);

xlabel('x/m')
ylabel('y/m')

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
    E=213e9;


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

    finray_info.E_A=LA/135e-3*E;
    finray_info.E_B=LB/156e-3*E;

    finray_info.constraint_ratio=constraint_ratio;

    finray_info.A_force_ratio=A_force_ratio;

    finray_info.B_force_ratio=B_force_ratio;


    FinrayA=finray_force(finray_info);




    x0=zeros(nA+nB+3+2*FinrayA.constraint_number,1);

    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

    f=@(x) FinrayA.cal_balance(x);
    [x_solve,fval,exitflag,output] = fsolve(f,x0,options);
    
    
    % 传感器的长度
    sensor_length=36e-3;
    
    sensor1_center=133.9e-3;
    sensor2_center=93.1e-3;


    p1=fix((sensor1_center-sensor_length/2-LB/nB/2)/(LB/nB))+1;
    q1=fix((sensor1_center+sensor_length/2-LB/nB/2)/(LB/nB))+1;

    p2=fix((sensor2_center-sensor_length/2-LB/nB/2)/(LB/nB))+1;
    q2=fix((sensor2_center+sensor_length/2-LB/nB/2)/(LB/nB))+1;



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
    plot(abs_pos(:,1),abs_pos(:,2),'o','MarkerSize',4,'LineWidth',2);
    axis equal
end