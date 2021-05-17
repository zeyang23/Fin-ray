% clear
% clc


% �޸���Լ�� r=30e-3
% ��1������
% delta=10e-3;
% Fx=-1.56;
% Fy=2.12;
% ��2������
% delta=15e-3;
% Fx=-2.20;
% Fy=3.42;
% ��3������
% delta=20e-3;
% Fx=-2.95;
% Fy=4.60;
% ��4������
% delta=25e-3;
% Fx=-3.00;
% Fy=5.80;
% ��5������
% delta=30e-3;
% Fx=-3.18;
% Fy=6.85;
% ��6������
% delta=35e-3;
% Fx=-3.40;
% Fy=7.82;
% ��7������
% delta=40e-3;
% Fx=-2.35;
% Fy=9.03;
% ��8������
% delta=45e-3;
% Fx=-3.73;
% Fy=9.82;
% ��9������
% delta=50e-3;
% Fx=-2.25;
% Fy=10.35;
% ��10������
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
