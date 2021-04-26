% ����cosserat model
% �������԰��shape sensing


clear
clc

L0=0.25;


% ����������
sensor_length=38e-3;

% 3����������λ��
S1=59e-3;
S2=125e-3;
S3=191e-3;


wid=14e-3;
thi=0.3e-3;
E=197*1e9;

Iz=1/12*thi.^3*wid;


load('coff_sensor_1.mat');
load('coff_sensor_2.mat');
load('coff_sensor_3.mat');

coff_sensor_1(2)=0;
coff_sensor_2(2)=0;
coff_sensor_3(2)=0;


U1_0=0;
U2_0=0;
U3_0=0;

U1=0;
U2=0;
U3=0;

U1_delta=U1-U1_0;
U2_delta=U2-U2_0;
U3_delta=U3-U3_0;

DeltaA=polyval(coff_sensor_1,U1_delta);
DeltaB=polyval(coff_sensor_2,U2_delta);
DeltaC=polyval(coff_sensor_3,U3_delta);


KA=DeltaA/sensor_length;
KB=DeltaB/sensor_length;
KC=DeltaC/sensor_length;


% �������ʷ������԰���״
g=@(x) check_shape(x,L0,E,Iz,KA,KB,KC,S1,S2,S3);


x0=zeros(6,1);

options_a = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
options_b = optimoptions('fsolve','Display','off','Algorithm','trust-region');
options_c = optimoptions('fsolve','Display','off','Algorithm','trust-region-dogleg');
[x_shape,fval_shape,exitflag_shape,output_shape]=fsolve(g,x0,options_a);

span = [0 L0];
y0 = [0;0;0;x_shape(1:3)];
options3=odeset('MaxStep',1e-2);
sol_Y2 = ode45(@(s,y) get_ydot(s,y,L0,E,Iz,x_shape(4:6)), span, y0,options3);

hold on
plot(sol_Y2.y(1,:),sol_Y2.y(2,:),'-b','MarkerSize',2)

axis equal


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

function res=check_balance(x,L,E,I,pdes)
    
    n0=x(1:2);
    m0=x(3);
    Fe=x(4:6);
    
    span = [0 L];
    y0 = [0;0;0;n0;m0];
    
    options=odeset('MaxStep',1e-2);
    
    [~,Y] = ode45(@(s,y) get_ydot(s,y,L,E,I,Fe), span, y0,options);
    
    ye=transpose(Y(end,:));
    
    res=zeros(6,1);
    res(1:2)=ye(1:2)-pdes(1:2);
    res(3)=ye(3)-pdes(3);
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