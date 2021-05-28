clear
clc

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

Delta=3:3:18;
hold on
title("接触力随夹钳移动的变化",'Color','blue','FontSize',14 )
xlabel('夹钳水平移动距离/mm')
ylabel('接触力/N')

plot(Delta,F_trans_series(3,:),'-o')
plot(Delta,Fcal_norm_series,'-*')
legend('实验测量','理论计算','Location','NorthWest')
axis([3 18 0 7])

print('force','-djpeg','-r1500')