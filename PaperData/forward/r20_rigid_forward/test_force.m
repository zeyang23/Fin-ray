clear
clc

Fsensor_series(:,1)=[0.36;1.45];
Fsensor_series(:,2)=[0.39;2.92];
Fsensor_series(:,3)=[0.43;4.65];
Fsensor_series(:,4)=[0.85;6.50];
Fsensor_series(:,5)=[1.03;9.36];
Fsensor_series(:,6)=[1.10;13.58];

Fcal_norm_series=[];
Fcal_norm_series(:,1)=1.51;
Fcal_norm_series(:,2)=2.95;
Fcal_norm_series(:,3)=4.54;
Fcal_norm_series(:,4)=6.44;
Fcal_norm_series(:,5)=9.10;
Fcal_norm_series(:,6)=13.58;

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

Delta=5:5:30;
hold on
title("接触力随夹钳移动的变化",'Color','blue','FontSize',14 )
xlabel('夹钳水平移动距离/mm')
ylabel('接触力/N')

plot(Delta,F_trans_series(3,:),'-o')
plot(Delta,Fcal_norm_series,'-*')
legend('实验测量','理论计算','Location','NorthWest')
axis([5 30 0 15])

print('force','-djpeg','-r1500')