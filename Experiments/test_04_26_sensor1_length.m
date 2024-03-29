% 弯曲传感器角度-电压标定
% 21-04-26

% 换算成传感器长度与电压的关系
% 注意长度的单位是mm

clear
clc

sensor_length=38;
H=12;

%%
volt1_1 = [];
volt1_1(1,:) = [100000000000000,2.02];
volt1_1(2,:) = [272,2.13];
volt1_1(3,:) = [203,2.17];
volt1_1(4,:) = [162,2.22];
volt1_1(5,:) = [134,2.28];
volt1_1(6,:) = [115,2.32];

volt1_1(7,:) = [100,2.38];
volt1_1(8,:) = [89,2.43];
volt1_1(9,:) = [80,2.49];
volt1_1(10,:) = [73,2.55];
volt1_1(11,:) = [67,2.61];

volt1_1(12,:) = [62,2.67];
volt1_1(13,:) = [57,2.74];
volt1_1(14,:) = [53,2.79];
volt1_1(15,:) = [50,2.85];
volt1_1(16,:) = [47,2.91];

% volt1_1(17,:) = [44,2.96];
% volt1_1(18,:) = [42,3.01];
% volt1_1(19,:) = [40,3.08];


% 减掉初值
init_1=volt1_1(1,2);
volt1_1(:,2)=volt1_1(:,2)-init_1;


for i=1:length(volt1_1)
    volt1_1(i,1)=get_length(sensor_length,volt1_1(i,1),H);
end
    
plot(volt1_1(:,2),volt1_1(:,1),'or','linewidth',1);

coff_1 = polyfit(volt1_1(:,2),volt1_1(:,1),1);

hold on
grid on
y_1 = polyval(coff_1,volt1_1(:,2));
plot(volt1_1(:,2),y_1,'--r','linewidth',1)


%%
volt1_2 = [];
volt1_2(1,:) = [100000000000000,2.09];
volt1_2(2,:) = [272,2.20];
volt1_2(3,:) = [203,2.24];
volt1_2(4,:) = [162,2.28];
volt1_2(5,:) = [134,2.33];
volt1_2(6,:) = [115,2.39];

volt1_2(7,:) = [100,2.45];
volt1_2(8,:) = [89,2.50];
volt1_2(9,:) = [80,2.56];
volt1_2(10,:) = [73,2.62];
volt1_2(11,:) = [67,2.67];

volt1_2(12,:) = [62,2.74];
volt1_2(13,:) = [57,2.80];
volt1_2(14,:) = [53,2.86];
volt1_2(15,:) = [50,2.91];
volt1_2(16,:) = [47,2.94];

% volt1_2(17,:) = [44,3.01];
% volt1_2(18,:) = [42,3.08];
% volt1_2(19,:) = [40,3.12];


% 减掉初值
init_2=volt1_2(1,2);
volt1_2(:,2)=volt1_2(:,2)-init_2;


for i=1:length(volt1_2)
    volt1_2(i,1)=get_length(sensor_length,volt1_2(i,1),H);
end

coff_2 = polyfit(volt1_2(:,2),volt1_2(:,1),1);

hold on
grid on
y_2 = polyval(coff_2,volt1_2(:,2));
plot(volt1_2(:,2),y_2,'--g','linewidth',1)


%%
coff_sensor1_length=1/2*(coff_1+coff_2);
save('coff_sensor1_length.mat','coff_sensor1_length');