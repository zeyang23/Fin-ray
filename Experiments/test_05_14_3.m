% 坐标变换


%% 柔性板坐标系变换到机械臂工具坐标系
clear
clc

theta=60/180*pi;


x_series=linspace(0,-100,11);
y_series=linspace(0,100,11);

p_tool_series=zeros(2,length(x_series));
for i=1:length(x_series)
    x=x_series(i);
    y=y_series(i);
    p_tool_series(:,i)=inv([sin(theta) cos(theta);cos(theta) -sin(theta)])*[x;y];
end




%% 力传感器坐标系变换到柔性板坐标系
clear
clc

theta=60/180*pi;

Fx_sensor=9.45;
Fy_sensor=-27.45;
Fx_trans=Fx_sensor*cos(theta)-Fy_sensor*sin(theta)
Fy_trans=-Fx_sensor*sin(theta)-Fy_sensor*cos(theta)
norm([Fx_sensor;Fy_sensor])