% Experiment-04-26
% Planar Rod Shape Sensing
% Algorithm: Cosserat

clear
clc

L0=0.25;


% 传感器长度
sensor_length=38e-3;

% 传感器中心的位置

% 从固定端到活动端，依次命名为A B C
SA=59e-3;
SB=125e-3;
SC=191e-3;


wid=14e-3;
thi=0.3e-3;
E=197*1e9;

Iz=1/12*thi.^3*wid;


load('coff_sensor_1.mat');
load('coff_sensor_2.mat');
load('coff_sensor_3.mat');

% coff_sensor_1(2)=0;
% coff_sensor_2(2)=0;
% coff_sensor_3(2)=0;


% 零位时的电压
U1_0=1.93;
U2_0=2.79;
U3_0=1.40;

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

% 第6组数据
U1=1.96;
U2=2.77;
U3=1.45;


U1_delta=U1-U1_0;
U2_delta=U2-U2_0;
U3_delta=U3-U3_0;

DeltaA=polyval(coff_sensor_3,U3_delta);
DeltaB=polyval(coff_sensor_1,U1_delta);
DeltaC=polyval(coff_sensor_2,U2_delta);


KA=DeltaA/sensor_length;
KB=DeltaB/sensor_length;
KC=DeltaC/sensor_length;


% 根据曲率反解柔性板形状
g=@(x) check_shape(x,L0,E,Iz,KA,KB,KC,SA,SB,SC);


x0=zeros(6,1);

options_a = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
options_b = optimoptions('fsolve','Display','off','Algorithm','trust-region');
options_c = optimoptions('fsolve','Display','off','Algorithm','trust-region-dogleg');
[x_shape,fval_shape,exitflag_shape,output_shape]=fsolve(g,x0,options_a);

span = [0 L0];
y0 = [0;0;0;x_shape(1:3)];
options3=odeset('MaxStep',1e-3);
sol_Y = ode45(@(s,y) get_ydot(s,y,L0,E,Iz,x_shape(4:6)), span, y0,options3);

hold on
plot(sol_Y.y(1,:),sol_Y.y(2,:),'-b','MarkerSize',2)

axis equal

YA=deval(sol_Y,SA);
YB=deval(sol_Y,SB);
YC=deval(sol_Y,SC);

plot(YA(1),YA(2),'o')
plot(YB(1),YB(2),'o')
plot(YC(1),YC(2),'o')

axis([0 L0 0 3/5*L0]);
axis equal

hold on

function res=check_shape(x,L,E,I,K1,K2,K3,S1,S2,S3)
    n0=x(1:2);
    m0=x(3);
    Fe=x(4:6);
    
    span = [0 L];
    y0 = [0;0;0;n0;m0];
    
    options=odeset('MaxStep',1e-2);
    
    sol_Y = ode45(@(s,y) get_ydot(s,y,L,E,I,Fe), span, y0,options);
    
    ye=sol_Y.y(:,end);
    
    Y1=deval(sol_Y,S1);
    Y2=deval(sol_Y,S2);
    Y3=deval(sol_Y,S3);

    
    res=zeros(6,1);
    res(1)=Y1(6)/E/I-K1;
    res(2)=Y2(6)/E/I-K2;
    res(3)=Y3(6)/E/I-K3;
    res(4:5)=ye(4:5)-Fe(1:2);
    res(6)=transpose(ye(1:2))*[0 1;-1 0]*ye(4:5)+ye(6)-transpose(ye(1:2))*[0 1;-1 0]*Fe(1:2)-Fe(3);

end


function ydot=get_ydot(s,y,L,E,I,Fe)
    ydot=zeros(size(y));
    
    r=y(1:2);
    theta=y(3);
    n=y(4:5);
    m=y(6);
    
%     delta=1e-16;
%     if (L-s)<delta
%         f=Fe(1:2);
%         l=Fe(3);
%     else
%         f=[0;0];
%         l=0;
%     end
    
    f=[0;0];
    l=0;

    ydot(1:2)=[cos(theta);sin(theta)];
    ydot(3)=1/(E*I)*m;
    ydot(4:5)=-f;
    ydot(6)=-[cos(theta) sin(theta)]*[0 1;-1 0]*n-l;
end

function [r,J]=cal_balance(x,Rod)
    Rod.theta=x(1:Rod.n_seg);
    Rod.F=x(Rod.n_seg+1:Rod.n_seg+3);
    
    Rod.update;
    r=Rod.cal_r;
    J=Rod.cal_rdot;
end