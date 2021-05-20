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
% 这组数据可以使用

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


% 有刚性约束 r=30e-3 31.5
% 第1组数据
% 这组数据可以使用
% Fsensor_series(:,1)=[0.30;1.15];
% Fsensor_series(:,2)=[0.37;2.46];
% Fsensor_series(:,3)=[0.30;4.19];
% Fsensor_series(:,4)=[0.25;6.34];
% Fsensor_series(:,5)=[0.35;9.00];
% Fsensor_series(:,6)=[0.80;12.60];

% Fcal_norm_series=[];
% Fcal_norm_series(:,1)=1.24;
% Fcal_norm_series(:,2)=2.67;
% Fcal_norm_series(:,3)=4.24;
% Fcal_norm_series(:,4)=6.06;
% Fcal_norm_series(:,5)=8.33;
% Fcal_norm_series(:,6)=12.81;



% 无刚性约束 r=30e-3 delta=40e-3
% 这组数据可以使用
% Fsensor_series(:,1)=[0.58;0.94];
% Fsensor_series(:,2)=[0.86;2.19];
% Fsensor_series(:,3)=[1.24;3.64];
% Fsensor_series(:,4)=[1.56;5.14];
% Fsensor_series(:,5)=[2.10;6.40];
% Fsensor_series(:,6)=[2.70;7.59];
% Fsensor_series(:,7)=[3.28;8.57];
% Fsensor_series(:,8)=[4.19;9.09];
% 
% Fcal_norm_series=[];
% Fcal_norm_series(:,1)=1.12;
% Fcal_norm_series(:,2)=2.57;
% Fcal_norm_series(:,3)=3.90;
% Fcal_norm_series(:,4)=5.54;
% Fcal_norm_series(:,5)=6.57;
% Fcal_norm_series(:,6)=7.57;
% Fcal_norm_series(:,7)=8.40;
% Fcal_norm_series(:,8)=9.08;
% 
% alpha=110.2/180*pi;


% 有刚性约束 椭圆a=20e-3 b=30e-3 21.5 31.5
% 第1组数据 0:5:25
% Fsensor_series(:,1)=[1.11;-0.28];
% Fsensor_series(:,2)=[2.58;-0.18];
% Fsensor_series(:,3)=[4.60;0.34];
% Fsensor_series(:,4)=[7.16;0.90];
% Fsensor_series(:,5)=[11.20;1.94];
% Fsensor_series(:,6)=[17.40;2.61];
% 
% 
% Fcal_norm_series=[];
% Fcal_norm_series(:,1)=1.23;
% Fcal_norm_series(:,2)=2.65;
% Fcal_norm_series(:,3)=3.90;
% Fcal_norm_series(:,4)=4.24;
% Fcal_norm_series(:,5)=6.50;
% Fcal_norm_series(:,6)=8.92;

% 第2组数据 0:3:18
% 这组数据可以使用
Fsensor_series(:,1)=[0.70;-0.24];
Fsensor_series(:,2)=[1.40;-0.42];
Fsensor_series(:,3)=[2.22;-0.42];
Fsensor_series(:,4)=[3.26;-0.21];
Fsensor_series(:,5)=[4.52;0.11];
Fsensor_series(:,6)=[6.04;0.52];


Fcal_norm_series=[];
Fcal_norm_series(:,1)=0.70;
Fcal_norm_series(:,2)=1.51;
Fcal_norm_series(:,3)=2.36;
Fcal_norm_series(:,4)=3.27;
Fcal_norm_series(:,5)=4.41;
Fcal_norm_series(:,6)=5.50;

alpha=200.2/180*pi;




F_trans_series=[];
for i=1:size(Fsensor_series,2)
    Fx=Fsensor_series(1,i);
    Fy=Fsensor_series(2,i);
    
    Fx_trans=-(Fx*cos(alpha)-Fy*sin(alpha));
    Fy_trans=-(Fx*sin(alpha)+Fy*cos(alpha));
    F_sensor_norm=norm([Fx;Fy]);
    
    F_trans_series(:,i)=[Fx_trans;Fy_trans;F_sensor_norm];
end