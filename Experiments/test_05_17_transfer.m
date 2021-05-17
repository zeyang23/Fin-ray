% clear
% clc


% 无刚性约束 r=30e-3
% 第1组数据
% delta=10e-3;
% Fx=-1.56;
% Fy=2.12;
% 第2组数据
% delta=15e-3;
% Fx=-2.20;
% Fy=3.42;
% 第3组数据
% delta=20e-3;
% Fx=-2.95;
% Fy=4.60;
% 第4组数据
% delta=25e-3;
% Fx=-3.00;
% Fy=5.80;
% 第5组数据
% delta=30e-3;
% Fx=-3.18;
% Fy=6.85;
% 第6组数据
% delta=35e-3;
% Fx=-3.40;
% Fy=7.82;
% 第7组数据
% delta=40e-3;
% Fx=-2.35;
% Fy=9.03;
% 第8组数据
% delta=45e-3;
% Fx=-3.73;
% Fy=9.82;
% 第9组数据
% delta=50e-3;
% Fx=-2.25;
% Fy=10.35;
% 第10组数据
% delta=55e-3;
% Fx=-3.40;
% Fy=11.12;

alpha=49.2/180*pi;

Fx_trans=-(Fx*cos(alpha)-Fy*sin(alpha));
Fy_trans=-(Fx*sin(alpha)+Fy*cos(alpha));

F_sensor_norm=norm([Fx;Fy]);

disp(Fx_trans)
disp(Fy_trans)
disp(F_sensor_norm)
