% clear
clc

Fsensor_series=[];

% r=30e-3

% Fsensor_series(:,1)=[-0.14;0.48];
% Fsensor_series(:,2)=[-0.20;1.25];
% Fsensor_series(:,3)=[-0.20;2.10];
% Fsensor_series(:,4)=[-0.90;2.93];
% Fsensor_series(:,5)=[-1.68;3.84];
% Fsensor_series(:,6)=[-2.55;4.96];
% Fsensor_series(:,7)=[-3.60;6.27];
% Fsensor_series(:,8)=[-4.44;7.80];
% Fsensor_series(:,9)=[-6.40;9.70];
% Fsensor_series(:,10)=[-7.94;11.96];

% alpha=69.8/180*pi;

% Fcal_series(:,1)=[0.707;1.12];
% Fcal_series(:,2)=[0.64;1.564];
% Fcal_series(:,3)=[0.6194;1.9169];
% Fcal_series(:,4)=[0.5873;2.4602];
% Fcal_series(:,5)=[0.5562;2.9587];
% Fcal_series(:,6)=[0.5535;3.6153];





% r=20e-3;

% Fsensor_series(:,1)=[-0.18;0.82];
% Fsensor_series(:,2)=[-0.26;1.61];
% Fsensor_series(:,3)=[-0.66;2.42];
% Fsensor_series(:,4)=[-1.25;3.28];
% Fsensor_series(:,5)=[-2.38;4.27];
% Fsensor_series(:,6)=[-3.42;5.50];
% Fsensor_series(:,7)=[-4.19;6.91];
% Fsensor_series(:,8)=[-5.22;8.62];
% Fsensor_series(:,9)=[-6.41;10.67];
% Fsensor_series(:,10)=[-6.00;12.40];
% 
% alpha=69.8/180*pi;
% 
% Fcal_series(:,1)=[0.6452;1.1904];
% Fcal_series(:,2)=[0.6014;1.6524];
% Fcal_series(:,3)=[0.5553;2.1058];
% Fcal_series(:,4)=[0.5574;2.7504];
% Fcal_series(:,5)=[0.5460;3.3737];
% Fcal_series(:,6)=[0.5435;4.1740];
% Fcal_series(:,7)=[0.5398;5.1303];
% Fcal_series(:,8)=[0.5500;6.2904];
% Fcal_series(:,9)=[0.5500;7.6037];
% Fcal_series(:,10)=[0.5492;8.8599];



% a20b30
Fsensor_series(:,1)=[-0.67;-0.20];
Fsensor_series(:,2)=[-1.40;-0.59];
Fsensor_series(:,3)=[-2.31;-1.06];
Fsensor_series(:,4)=[-3.24;-2.05];
Fsensor_series(:,5)=[-4.45;-3.70];
Fsensor_series(:,6)=[-5.79;-4.85];
Fsensor_series(:,7)=[-7.58;-6.98];
Fsensor_series(:,8)=[-9.59;-8.84];
alpha=-20.2;

Fcal_series(:,1)=[0.6755;1.2464];
Fcal_series(:,2)=[0.6338;1.6597];
Fcal_series(:,3)=[0.6117;2.1036];
Fcal_series(:,4)=[0.5855;2.5608];
Fcal_series(:,5)=[0.5551;3.1720];
Fcal_series(:,6)=[0.5573;3.9010];


F_trans_series=[];
for i=1:size(Fsensor_series,2)
    Fx=Fsensor_series(1,i);
    Fy=Fsensor_series(2,i);
    
    Fx_trans=-(Fx*cos(alpha)-Fy*sin(alpha));
    Fy_trans=-(Fx*sin(alpha)+Fy*cos(alpha));
    F_sensor_norm=norm([Fx;Fy]);
    
    F_trans_series(:,i)=[Fx_trans;Fy_trans;F_sensor_norm];
end