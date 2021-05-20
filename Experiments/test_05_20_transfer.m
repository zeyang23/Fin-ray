% clear
clc

Fsensor_series=[];

% 有刚性约束 r=20e-3
% delta 5~30
% 第1组数据
% Fsensor_series(:,1)=[0.35;1.00];
% Fsensor_series(:,2)=[0.40;2.25];
% Fsensor_series(:,3)=[-0.04;4.07];
% Fsensor_series(:,4)=[-0.18;6.14];
% Fsensor_series(:,5)=[-0.20;8.76];
% Fsensor_series(:,6)=[-0.05;11.76];


% 第2组数据
% Fsensor_series(:,1)=[0.11;0.38];
% Fsensor_series(:,2)=[0.47;1.55];
% Fsensor_series(:,3)=[0.35;3.08];
% Fsensor_series(:,4)=[0.42;4.82];
% Fsensor_series(:,5)=[0.40;7.10];
% Fsensor_series(:,6)=[0.60;10.12];
% Fsensor_series(:,7)=[1.22;13.75];

% 第3组数据 这组数据比较好
% 错位对齐，1mm的误差可以有很大的影响

% Fsensor_series(:,1)=[0.36;1.45];
% Fsensor_series(:,2)=[0.39;2.92];
% Fsensor_series(:,3)=[0.43;4.65];
% Fsensor_series(:,4)=[0.85;6.50];
% Fsensor_series(:,5)=[1.03;9.36];
% Fsensor_series(:,6)=[1.10;13.58];
% 
% Fcal_norm_series=[];
% Fcal_norm_series(:,1)=1.51;
% Fcal_norm_series(:,2)=2.95;
% Fcal_norm_series(:,3)=4.54;
% Fcal_norm_series(:,4)=6.44;
% Fcal_norm_series(:,5)=9.10;
% Fcal_norm_series(:,6)=13.58;


% 有刚性约束 r=30e-3
% 第1组数据
Fsensor_series(:,1)=[0.30;1.15];
Fsensor_series(:,2)=[0.37;2.46];
Fsensor_series(:,3)=[0.30;4.19];
Fsensor_series(:,4)=[0.25;6.34];
Fsensor_series(:,5)=[0.35;9.00];
Fsensor_series(:,6)=[0.80;12.60];

% Fcal_norm_series=[];
Fcal_norm_series(:,1)=1.24;
Fcal_norm_series(:,2)=2.67;
Fcal_norm_series(:,3)=4.24;
Fcal_norm_series(:,4)=6.06;
Fcal_norm_series(:,5)=8.33;
Fcal_norm_series(:,6)=12.81;



alpha=110.2/180*pi;

F_trans_series=[];
for i=1:size(Fsensor_series,2)
    Fx=Fsensor_series(1,i);
    Fy=Fsensor_series(2,i);
    
    Fx_trans=-(Fx*cos(alpha)-Fy*sin(alpha));
    Fy_trans=-(Fx*sin(alpha)+Fy*cos(alpha));
    F_sensor_norm=norm([Fx;Fy]);
    
    F_trans_series(:,i)=[Fx_trans;Fy_trans;F_sensor_norm];
end