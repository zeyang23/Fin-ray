% Experiment-04-26
% Planar Rod Shape Sensing
% Algorithm: Principle Axis Decomposition

clear
clc

wid=14e-3;
thi=0.3e-3;
E=197*1e9;

L0=0.25;
n=100;

% 传感器的长度
sensor_length=38e-3;
sensor_length_ratio=sensor_length/L0;

% 传感器的中心位置

% 从固定端到活动端，依次命名为A B C

sensorA_center_ratio=59e-3/L0;
sensorB_center_ratio=125e-3/L0;
sensorC_center_ratio=191e-3/L0;


pA=fix((sensorA_center_ratio-sensor_length_ratio/2)*n);
qA=fix((sensorA_center_ratio+sensor_length_ratio/2)*n);

pB=fix((sensorB_center_ratio-sensor_length_ratio/2)*n);
qB=fix((sensorB_center_ratio+sensor_length_ratio/2)*n);

pC=fix((sensorC_center_ratio-sensor_length_ratio/2)*n);
qC=fix((sensorC_center_ratio+sensor_length_ratio/2)*n);



DA=zeros(1,n);
for i=pA:qA
    DA(i)=1;
end

DB=zeros(1,n);
for i=pB:qB
    DB(i)=1;
end

DC=zeros(1,n);
for i=pC:qC
    DC(i)=1;
end


%% 使用fsolve求解

x0=zeros(n+3,1);
options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');


load('coff_sensor_1.mat');
load('coff_sensor_2.mat');
load('coff_sensor_3.mat');

% coff_sensor_1(2)=0;
% coff_sensor_2(2)=0;
% coff_sensor_3(2)=0;


%% 第1次实验
% % 零位时的电压
% U1_0=1.93;
% U2_0=2.79;
% U3_0=1.40;

% % 第1组数据
% U1=2.01;
% U2=3.00;
% U3=1.43;

% % 第2组数据
% U1=2.12;
% U2=3.15;
% U3=1.45;

% % 第3组数据
% U1=2.11;
% U2=3.30;
% U3=1.44;

% % 第4组数据
% U1=1.97;
% U2=2.95;
% U3=1.45;

% % 第5组数据
% U1=2.00;
% U2=2.79;
% U3=1.50;

% % 第6组数据
% U1=1.96;
% U2=2.77;
% U3=1.45;


%% 第2次实验
% % 零位时的电压
% U1_0=1.97;
% U2_0=2.78;
% U3_0=1.40;

% % 第1组数据
% U1=2.05;
% U2=2.84;
% U3=1.43;

% % 第2组数据
% U1=2.04;
% U2=2.94;
% U3=1.40;

% % 第3组数据
% U1=2.05;
% U2=2.98;
% U3=1.43;

% % 第4组数据
% U1=2.08;
% U2=3.02;
% U3=1.47;

% % 第5组数据
% U1=2.13;
% U2=3.18;
% U3=1.55;

% % 第6组数据
% U1=2.09;
% U2=2.83;
% U3=1.47;

% % 第7组数据
% U1=1.99;
% U2=2.79;
% U3=1.56;



%% 第3次实验
% % 零位时的电压
% U1_0=1.96;
% U2_0=2.78;
% U3_0=1.34;

% % 第1组数据
% U1=2.01;
% U2=2.83;
% U3=1.35;

% % 第2组数据
% U1=2.03;
% U2=2.83;
% U3=1.37;

% % 第3组数据
% U1=2.05;
% U2=2.95;
% U3=1.37;

% % 第4组数据
% U1=2.17;
% U2=3.18;
% U3=1.34;

% % 第5组数据
% U1=2.20;
% U2=3.09;
% U3=1.35;

% % 第6组数据
% U1=2.13;
% U2=3.02;
% U3=1.36;

% % 第7组数据
% U1=2.35;
% U2=3.30;
% U3=1.35;

% % 第8组数据
% U1=2.25;
% U2=3.55;
% U3=1.35;


%% 第4次实验
% % 零位时的电压
% U1_0=1.97;
% U2_0=2.78;
% U3_0=1.40;

% % 第1组数据
% U1=2.47;
% U2=3.51;
% U3=1.42;

% % 第2组数据
% U1=2.29;
% U2=3.11;
% U3=1.42;


% % 第3组数据
% U1=2.59;
% U2=3.57;
% U3=1.50;


% % 第4组数据
% U1=2.97;
% U2=3.86;
% U3=1.54;


%%

U1_delta=U1-U1_0;
U2_delta=U2-U2_0;
U3_delta=U3-U3_0;

DeltaA=polyval(coff_sensor_3,U3_delta);
DeltaB=polyval(coff_sensor_1,U1_delta);
DeltaC=polyval(coff_sensor_2,U2_delta);


f_sensor=@(x) myfun(x,E,L0,wid,thi,n,DA,DB,DC,DeltaA,DeltaB,DeltaC);
[x_sensor,fval_sensor,exitflag_sensor,output_sensor] = fsolve(f_sensor,x0,options);


Rod_solve=planar_nR(E,L0,wid,thi,n,[0;0;0]);
Rod_solve.theta=x_sensor(1:n);
Rod_solve.F=x_sensor(end-2:end);

Rod_solve.update;
Rod_solve.plot_all;

Rod_solve.cal_posall;
sensor1_pos=Rod_solve.pos_all(1+pA:1+qA,:);
sensor2_pos=Rod_solve.pos_all(1+pB:1+qB,:);
sensor3_pos=Rod_solve.pos_all(1+pC:1+qC,:);

hold on
plot(sensor1_pos(:,1),sensor1_pos(:,2),'o')
plot(sensor2_pos(:,1),sensor2_pos(:,2),'o')
plot(sensor3_pos(:,1),sensor3_pos(:,2),'o')

axis([0 L0 0 3/5*L0]);
axis equal

hold on

Delta1_degree=DeltaA/pi*180;
Delta2_degree=DeltaB/pi*180;
Delta3_degree=DeltaC/pi*180;

disp("第1个夹角")
disp(Delta1_degree)
disp("第2个夹角")
disp(Delta2_degree)
disp("第3个夹角")
disp(Delta3_degree)


function [f,grad] = myfun(x,E,L0,wid,thi,n,D1,D2,D3,Delta1,Delta2,Delta3)
    Rod=planar_nR(E,L0,wid,thi,n,[0;0;0]);
    
    theta=zeros(n,1);
    F=zeros(3,1);
    
    theta=x(1:n);
    F=x(end-2:end);
    
    Rod.theta=theta;
    Rod.F=F;
    
    Rod.update;
    
    r1=Rod.K_theta*theta-transpose(Rod.Jacobian)*F;
    r2=D1*theta-Delta1;
    r3=D2*theta-Delta2;
    r4=D3*theta-Delta3;
    
    f=[r1;r2;r3;r4];
    
    if nargout > 1
        J=zeros(n+3);
    
        J(1:n,1:n)=Rod.K_theta-1*Rod.partial;
        J(1:n,n+1:end)=-transpose(Rod.Jacobian);
        J(n+1,1:n)=D1;
        J(n+2,1:n)=D2;
        J(n+3,1:n)=D3;
        grad=J;
    end
end