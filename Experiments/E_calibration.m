% 使用力传感器测定末端力，测定多组数据后使用最小二乘法标定杨氏模量

% 杨氏模量的单位需要取Mpa
% 如果取为pa，则绝对数值过大，fmincon无法求解
% 使用lsqnonlin比使用fmincon的效果更好

clear
clc

Rod_info.wid=14;
Rod_info.thi=0.2;

L0=350;
Rod_info.L0=L0;
Rod_info.n=100;


pdes_series=[];
pdes_series(:,1)=[L0-0.4*L0;0.4*L0;pi/3];
pdes_series(:,2)=[L0-0.3*L0;0.3*L0;pi/6];
pdes_series(:,3)=[L0-0.2*L0;0.2*L0;pi/12];

Fsensor_series=[];
Fsensor_series(:,1)=[-0.3887;-0.2507;-0.0245];
Fsensor_series(:,2)=[-0.5136;-0.2104;-0.0291];
Fsensor_series(:,3)=[-0.5670;-0.1376;-0.0264];


calibration_resnorm=@(E) cal_resnorm(E,Rod_info,pdes_series,Fsensor_series);
calibration_res=@(E) cal_res(E,Rod_info,pdes_series,Fsensor_series);

E0=190;

% [E_calib_fmincon,resnorm_fmincon]=fmincon(calibration_resnorm,E0);
[E_calib_lsqnonlin,resnorm_lsqnonlin]=lsqnonlin(calibration_res,E0);





function resnorm=cal_resnorm(E,Rod_info,pdes_series,Fsensor_series)
    number=size(pdes_series,2);
    res=zeros(3*number);
    
    for i=1:number
        Rod=planar_nR(E,Rod_info.L0,Rod_info.wid,Rod_info.thi,Rod_info.n,pdes_series(:,i));
        
        f=@(x) cal_balance(x,Rod);
        options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');
        x0=zeros(Rod_info.n+3,1);
        [x_solve,~,exitflag,~] = fsolve(f,x0,options);
        if exitflag<=0
            error('fail')
        end
        res(3*i-2:3*i-1)=x_solve(Rod.n_seg+1:Rod.n_seg+2)*1000-Fsensor_series(1:2,i);
    end
    
    resnorm=norm(res);
        
end

function res=cal_res(E,Rod_info,pdes_series,Fsensor_series)
    number=size(pdes_series,2);
    res=zeros(3*number);
    
    for i=1:number
        Rod=planar_nR(E,Rod_info.L0,Rod_info.wid,Rod_info.thi,Rod_info.n,pdes_series(:,i));
        
        f=@(x) cal_balance(x,Rod);
        options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');
        x0=zeros(Rod_info.n+3,1);
        [x_solve,~,exitflag,~] = fsolve(f,x0,options);
        if exitflag<=0
            error('fail')
        end
        res(3*i-2:3*i-1)=x_solve(Rod.n_seg+1:Rod.n_seg+2)*1000-Fsensor_series(1:2,i);
    end
        
end

function [r,J]=cal_balance(x,Rod)
    Rod.theta=x(1:Rod.n_seg);
    Rod.F=x(Rod.n_seg+1:Rod.n_seg+3);
    
    Rod.update;
    r=Rod.cal_r;
    J=Rod.cal_rdot;
end