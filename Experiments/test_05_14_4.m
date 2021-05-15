% 使用力传感器测定末端力，测定多组数据后使用最小二乘法标定杨氏模量

% 杨氏模量的单位需要取Gpa
% 如果取为pa，则绝对数值过大，无法求解
% 使用lsqnonlin比使用fmincon的效果更好

clear
clc

Rod_info.wid=14;
Rod_info.thi=0.5;

L0=250;
Rod_info.L0=L0;
Rod_info.n=50;


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

Fsensor_series=[];

% 第一次实验，对零位上的力做清零
% 经过校验，得出结论：不能瞎清零
% Fsensor_series(:,1)=[14.0;-26.85;0];
% Fsensor_series(:,2)=[13.55;-27.8;0];
% Fsensor_series(:,3)=[12.87;-28.65;0];
% Fsensor_series(:,4)=[12.04;-29.40;0];
% Fsensor_series(:,5)=[11.08;-30.10;0];
% Fsensor_series(:,6)=[9.98;-30.74;0];
% Fsensor_series(:,7)=[8.75;-31.28;0];
% Fsensor_series(:,8)=[7.43;-31.7;0];
% Fsensor_series(:,9)=[5.96;-32;0];
% Fsensor_series(:,10)=[4.41;-32.15;0];


% 第二次实验，对零位上的力不清零
% 经过校验，数据不错
init_zero=[-5.3;9.35;0];

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

final_zero=[-7.06;12.55;0];


calibration_res=@(E) cal_res(E,Rod_info,pdes_series,Fsensor_series);
calibration_Fnorm_res=@(E) cal_Fnorm_res(E,Rod_info,pdes_series,Fsensor_series);

E0=197;

[E_calib_lsqnonlin,resnorm_lsqnonlin,residual_lsqnonlin,~,~]=lsqnonlin(calibration_res,E0);
% [E_calib_lsqnonlin,resnorm_lsqnonlin,residual_lsqnonlin,~,~]=lsqnonlin(calibration_Fnorm_res,E0);


function res=cal_res(E,Rod_info,pdes_series,Fsensor_series)
    number=size(pdes_series,2);
    res=zeros(2*number,1);
    
    x0=zeros(Rod_info.n+3,1);
    
    for i=1:number
        Rod=planar_nR(E,Rod_info.L0,Rod_info.wid,Rod_info.thi,Rod_info.n,pdes_series(:,i));
        
        f=@(x) cal_balance(x,Rod);
        options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');
        [x_solve,~,exitflag,~] = fsolve(f,x0,options);
        if exitflag<=0
            error('fail')
        end
        
        Fx_sensor=Fsensor_series(1,i);
        Fy_sensor=Fsensor_series(2,i);
        theta=60/180*pi;
        Fx_sensor_transfered=Fx_sensor*cos(theta)-Fy_sensor*sin(theta);
        Fy_sensor_transfered=-Fx_sensor*sin(theta)-Fy_sensor*cos(theta);
        
        x0=x_solve;
        res(2*i-1:2*i)=x_solve(Rod.n_seg+1:Rod.n_seg+2)*1000+[Fx_sensor_transfered;Fy_sensor_transfered];
    end
        
end

function res=cal_Fnorm_res(E,Rod_info,pdes_series,Fsensor_series)
    x0=zeros(Rod_info.n+3,1);
    number=size(pdes_series,2);
    res=zeros(number,1);
    
    for i=1:number
        Rod=planar_nR(E,Rod_info.L0,Rod_info.wid,Rod_info.thi,Rod_info.n,pdes_series(:,i));
        
        f=@(x) cal_balance(x,Rod);
        options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');
        
        [x_solve,~,exitflag,~] = fsolve(f,x0,options);
        if exitflag<=0
            error('fail')
        end
        res(i)=norm(x_solve(Rod.n_seg+1:Rod.n_seg+2))*1000-norm(Fsensor_series(1:2,i));
        
        x0=x_solve;
    end
        
end

function [r,J]=cal_balance(x,Rod)
    Rod.theta=x(1:Rod.n_seg);
    Rod.F=x(Rod.n_seg+1:Rod.n_seg+3);
    
    Rod.update;
    r=Rod.cal_r;
    J=Rod.cal_rdot;
end