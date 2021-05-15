% 对比实验数据和理论解
clear
clc

wid=14;
thi=0.5;
E=211;

L0=255; % 柔性板的长度该怎么取？ 真正的长度应该在250和255之间
n=50;


pdes_series=[];
pdes_series(:,1)=[L0-10;10;0];
pdes_series(:,2)=[L0-20;20;0];
pdes_series(:,3)=[L0-30;30;0];
pdes_series(:,4)=[L0-40;40;0];
pdes_series(:,5)=[L0-50;50;0];
pdes_series(:,6)=[L0-60;60;0];
pdes_series(:,7)=[L0-70;70;0];
pdes_series(:,8)=[L0-80;80;0];
pdes_series(:,9)=[L0-90;90;0];
pdes_series(:,10)=[L0-100;100;0];

Fsensor_series(:,1)=[8.34;-16.7;0];
Fsensor_series(:,2)=[7.85;-17.68;0];
Fsensor_series(:,3)=[7.16;-18.5;0];
Fsensor_series(:,4)=[6.34;-19.3;0];
Fsensor_series(:,5)=[5.40;-20;0];
Fsensor_series(:,6)=[4.30;-20.63;0];
Fsensor_series(:,7)=[3.08;-21.20;0];
Fsensor_series(:,8)=[1.75;-21.61;0];
Fsensor_series(:,9)=[0.30;-21.92;0];
Fsensor_series(:,10)=[-1.24;-22.07;0];


x0=zeros(n+3,1);


F_theory=zeros(4,length(pdes_series));


for i=1:length(pdes_series)
    pdes=pdes_series(:,i);

    Rod=planar_nR(E,L0,wid,thi,n,pdes);


    f=@(x) cal_balance(x,Rod);
    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

    [x_solve,fval,exitflag,output] = fsolve(f,x0,options);
    
    F_theory(1:2,i)=1000*x_solve(Rod.n_seg+1:Rod.n_seg+2);
    F_theory(3,i)=norm(F_theory(1:2,i));
    F_theory(4,i)=atan2(F_theory(2,i),F_theory(1,i))/pi*180;
    
    x0=x_solve;
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

DELTA=F_sensor_transfered-F_theory;

function [r,J]=cal_balance(x,Rod)
    Rod.theta=x(1:Rod.n_seg);
    Rod.F=x(Rod.n_seg+1:Rod.n_seg+3);
    
    Rod.update;
    r=Rod.cal_r;
    J=Rod.cal_rdot;
end