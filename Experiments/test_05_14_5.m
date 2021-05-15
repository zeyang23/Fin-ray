% 第3次实验 从[0;0]到[-150;150]

% 对比实验数据和两种理论解

clear
clc

wid=14e-3;
thi=0.5e-3;

Iz=1/12*thi.^3*wid;

E=211e9; % 65Mn的杨氏模量到底是多少？ 查到的说法不一

L0=255e-3; % 柔性板的长度该怎么取？ 真正的长度应该在250和255之间
n=50;


pdes_series=[];
pdes_series(:,1)=[L0-10e-3;10e-3;0];
pdes_series(:,2)=[L0-20e-3;20e-3;0];
pdes_series(:,3)=[L0-30e-3;30e-3;0];
pdes_series(:,4)=[L0-40e-3;40e-3;0];
pdes_series(:,5)=[L0-50e-3;50e-3;0];
pdes_series(:,6)=[L0-60e-3;60e-3;0];
pdes_series(:,7)=[L0-70e-3;70e-3;0];
pdes_series(:,8)=[L0-80e-3;80e-3;0];
pdes_series(:,9)=[L0-90e-3;90e-3;0];
pdes_series(:,10)=[L0-100e-3;100e-3;0];
pdes_series(:,11)=[L0-110e-3;110e-3;0];
pdes_series(:,12)=[L0-120e-3;120e-3;0];
pdes_series(:,13)=[L0-130e-3;130e-3;0];
pdes_series(:,14)=[L0-140e-3;140e-3;0];
pdes_series(:,15)=[L0-150e-3;150e-3;0];

Fsensor_series(:,1)=[8.34;-16.8;0];
Fsensor_series(:,2)=[7.85;-17.78;0];
Fsensor_series(:,3)=[7.18;-18.59;0];
Fsensor_series(:,4)=[6.35;-19.35;0];
Fsensor_series(:,5)=[5.43;-20.06;0];
Fsensor_series(:,6)=[4.31;-20.70;0];
Fsensor_series(:,7)=[3.08;-21.22;0];
Fsensor_series(:,8)=[1.75;-21.65;0];
Fsensor_series(:,9)=[0.30;-21.95;0];
Fsensor_series(:,10)=[-1.25;-22.07;0];
Fsensor_series(:,11)=[-2.83;-22.04;0];
Fsensor_series(:,12)=[-4.45;-21.81;0];
Fsensor_series(:,13)=[-6.01;-21.37;0];
Fsensor_series(:,14)=[-7.50;-20.75;0];
Fsensor_series(:,15)=[-8.79;-19.94;0];


x0_pad=zeros(n+3,1);


F_theory_pad=zeros(4,length(pdes_series));

Fe_series=[];


for i=1:length(pdes_series)
    pdes=pdes_series(:,i);

    Rod=planar_nR(E,L0,wid,thi,n,pdes);


    f=@(x) cal_balance_pad(x,Rod);
    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

    [x_pad,fval,exitflag,output] = fsolve(f,x0_pad,options);
    
    F_theory_pad(1:2,i)=x_pad(Rod.n_seg+1:Rod.n_seg+2);
    F_theory_pad(3,i)=norm(F_theory_pad(1:2,i));
    F_theory_pad(4,i)=atan2(F_theory_pad(2,i),F_theory_pad(1,i))/pi*180;
    
    x0_pad=x_pad;
    
    Rod.theta=x_pad(1:Rod.n_seg);
    Rod.F=x_pad(Rod.n_seg+1:Rod.n_seg+3);
    
    Fe_series(:,i)=Rod.F;

    Rod.update;
    Rod.plot_all;
    hold on
end



F_theory_cosserat=zeros(4,length(pdes_series));

% cosserat的多解现象很严重，初值是个玄学
x0_cosserat=zeros(6,1);
x0_cosserat(1:3)=Fe_series(:,end);
x0_cosserat(4:6)=Fe_series(:,end);

for i=1:length(pdes_series)
    pdes=pdes_series(:,length(pdes_series)+1-i);
    f=@(x) check_balance_cosserat(x,L0,E,Iz,pdes);
    
%     x0_cosserat(1:3)=-Fe_series(:,length(pdes_series)+1-i)-[0;0;transpose(pdes(1:2))*[0 1;-1 0]*Fe_series(1:2,length(pdes_series)+1-i)];
%     x0_cosserat(4:6)=Fe_series(:,length(pdes_series)+1-i);

    options_a = optimoptions('fsolve','Display','off','Algorithm','trust-region');
    [x_cosserat,fval_cosserat,exitflag_cosserat,output_cosserat]=fsolve(f,x0_cosserat,options_a);
    
    F_theory_cosserat(1:2,length(pdes_series)+1-i)=x_cosserat(4:5);
    F_theory_cosserat(3,length(pdes_series)+1-i)=norm(F_theory_cosserat(1:2,length(pdes_series)+1-i));
    F_theory_cosserat(4,length(pdes_series)+1-i)=atan2(F_theory_cosserat(2,length(pdes_series)+1-i),F_theory_cosserat(1,length(pdes_series)+1-i))/pi*180;
    
    span = [0 L0];
    y0 = [0;0;0;x_cosserat(1:3)];
    options3=odeset('MaxStep',1e-3);
    [s,Y] = ode45(@(s,y) get_ydot(s,y,L0,E,Iz,x_cosserat(4:6)), span, y0,options3);

    hold on
    plot(Y(:,1),Y(:,2))

    x0_cosserat=x_cosserat;
end

F_sensor_transfered=zeros(4,length(Fsensor_series));
for i=1:length(Fsensor_series)
    Fx_sensor=Fsensor_series(1,i);
    Fy_sensor=Fsensor_series(2,i);
    theta=60/180*pi;
    F_sensor_transfered(1,i)=Fx_sensor*cos(theta)-Fy_sensor*sin(theta);
    F_sensor_transfered(2,i)=-Fx_sensor*sin(theta)-Fy_sensor*cos(theta);
    
    F_sensor_transfered(1:2,i)=-F_sensor_transfered(1:2,i);
    
    F_sensor_transfered(3,i)=norm(F_sensor_transfered(1:2,i));
    F_sensor_transfered(4,i)=atan2(F_sensor_transfered(2,i),F_sensor_transfered(1,i))/pi*180;
end

DELTA_pad=F_sensor_transfered-F_theory_pad;
DELTA_cosserat=F_sensor_transfered-F_theory_cosserat;

DELTA_pad_cosserat=F_theory_pad-F_theory_cosserat;


function [r,J]=cal_balance_pad(x,Rod)
    Rod.theta=x(1:Rod.n_seg);
    Rod.F=x(Rod.n_seg+1:Rod.n_seg+3);
    
    Rod.update;
    r=Rod.cal_r;
    J=Rod.cal_rdot;
end

function res=check_balance_cosserat(x,L,E,I,pdes)
    
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