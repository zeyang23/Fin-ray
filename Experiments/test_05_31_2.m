% 假设接触点的位置已知，根据背部曲率求解柔性夹钳形状
% 不使用shooting method，直接求解整个方程组

% 21-05-31
% 通过了梯度检测
% 对初值敏感
% 对算法敏感
% 例如 trust-region-dogleg就不太行



clear
clc

% r=30e-3

%% 根据传感器示数，求出夹角

coff_1=[1.186256198447459,0.056997358385406];
coff_2=[1.151462241145653,0.027617725599564];


U1_0=3.19-1.51;
U2_0=1.51;
% 第1组数据 delta=3
% U1=3.21-1.52;
% U2=1.52;

% 第2组数据 delta=6
% U1=3.27-1.55;
% U2=1.55;

% 第3组数据 delta=9
% U1=3.31-1.57;
% U2=1.57;

% 第4组数据 delta=12
U1=3.37-1.605;
U2=1.605;

% 第5组数据 delta=15
% U1=3.42-1.64;
% U2=1.64;

% 第6组数据 delta=18
% U1=3.48-1.67;
% U2=1.67;

% 第7组数据 delta=21
% U1=3.55-1.71;
% U2=1.71;

% % 第8组数据 delta=24
% U1=3.63-1.75;
% U2=1.75;

% 第9组数据 delta=27
% U1=3.70-1.79;
% U2=1.79;

% 第10组数据 delta=30
% U1=3.78-1.83;
% U2=1.83;



U1_delta=U1-U1_0;
U2_delta=U2-U2_0;


DeltaA=polyval(coff_1,U1_delta);
DeltaB=polyval(coff_2,U2_delta);


lambda1=0.45;
lambda2=0.45;



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
A_force_ratio=[0,lambda1;
               0,lambda2];


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


% 传感器的长度
sensor_length=36e-3;

% 传感器的中心位置
sensor1_center=133.9e-3;
sensor2_center=93.1e-3;


p1=fix((sensor1_center-sensor_length/2-LB/nB/2)/(LB/nB))+1;
q1=fix((sensor1_center+sensor_length/2-LB/nB/2)/(LB/nB))+1;

p2=fix((sensor2_center-sensor_length/2-LB/nB/2)/(LB/nB))+1;
q2=fix((sensor2_center+sensor_length/2-LB/nB/2)/(LB/nB))+1;




x0=zeros(nA+nB+3+2*FinrayA.constraint_number+2,1);

f=@(x) cal_all_balance(x,FinrayA,p1,q1,p2,q2,DeltaA,DeltaB);

options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off','Algorithm','levenberg-marquardt');

[x_solve,fval,exitflag,output] = fsolve(f,x0,options);

F_solve=sum(x_solve(end-1:end));
resnorm=norm(fval)


%%

figure
FinrayA.plot_state(x_solve);


FinrayA.RodB.cal_posall;
sensor1_pos=FinrayA.RodB.pos_all(1+p1:1+q1,:);
sensor2_pos=FinrayA.RodB.pos_all(1+p2:1+q2,:);

plot_sensor_pos(sensor1_pos,beta,[xB,yB]);
plot_sensor_pos(sensor2_pos,beta,[xB,yB]);


axis equal
axis([-80e-3 100e-3 0 160e-3]);

xlabel('x/m')
ylabel('y/m')

rigid_abs_pos=FinrayA.get_rigid_pos(x_solve);
rigid_abs_pos(:,1)=rigid_abs_pos(:,1)-xA;
rigid_abs_pos=1000*rigid_abs_pos;

% print('r30_delta18','-djpeg','-r600')
% save('r30_delta18_cal.mat','rigid_abs_pos')


function [r,J]=cal_all_balance(x,Finray,p1,q1,p2,q2,Delta1,Delta2)

    Finray.A_force_index(1,1)=x(end-1);
    Finray.A_force_index(2,1)=x(end);
    
    Finray.A_force_ratio(1,1)=x(end-1);
    Finray.A_force_ratio(2,1)=x(end);
    
    [r1,J1]=Finray.cal_balance(x(1:end-2));
    
    r=zeros(length(r1)+2,1);
    J=zeros(length(r1)+2);
    
    r(1:length(r1))=r1;
    J(1:length(r1),1:length(r1))=J1;
    
    Delta1_now=sum(Finray.RodB.theta(p1:q1));
    Delta2_now=sum(Finray.RodB.theta(p2:q2));
    
    r(length(r1)+1)=Delta1_now-Delta1;
    r(length(r1)+2)=Delta2_now-Delta2;
    
    D1=zeros(1,Finray.RodB.n_seg);
    for i=p1:q1
        D1(i)=1;
    end
    
    D2=zeros(1,Finray.RodB.n_seg);
    for i=p2:q2
        D2(i)=1;
    end
    
    
    
    JF1=Finray.A_force_array(1).Jacobian;
    JF2=Finray.A_force_array(2).Jacobian;
    
    lambdaF1=zeros(Finray.A_force_index(1,2),Finray.nA);
    lambdaF1(1:Finray.A_force_index(1,2),1:Finray.A_force_index(1,2))=diag(ones(Finray.A_force_index(1,2),1));
    
    lambdaF2=zeros(Finray.A_force_index(2,2),Finray.nA);
    lambdaF2(1:Finray.A_force_index(2,2),1:Finray.A_force_index(2,2))=diag(ones(Finray.A_force_index(2,2),1));
    
    pk_F1=Finray.A_force_array(1).pe;
    pk_F2=Finray.A_force_array(2).pe;
    
    J(4:Finray.RodA.n_seg+3,end-1)=-transpose(lambdaF1)*transpose(JF1)*[sin(pk_F1(3));-cos(pk_F1(3));0];
    J(4:Finray.RodA.n_seg+3,end)=-transpose(lambdaF2)*transpose(JF2)*[sin(pk_F2(3));-cos(pk_F2(3));0];
    
    J(length(r1)+1,Finray.nA+1:Finray.nA+Finray.nB)=D1;
    J(length(r1)+2,Finray.nA+1:Finray.nA+Finray.nB)=D2;
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