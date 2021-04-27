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
volt1_1(1,:) = [100000000000000,1.40];
volt1_1(2,:) = [272,1.49];
volt1_1(3,:) = [203,1.53];
volt1_1(4,:) = [162,1.57];
volt1_1(5,:) = [134,1.61];
volt1_1(6,:) = [115,1.65];

volt1_1(7,:) = [100,1.70];
volt1_1(8,:) = [89,1.74];
volt1_1(9,:) = [80,1.78];
volt1_1(10,:) = [73,1.82];
volt1_1(11,:) = [67,1.88];

volt1_1(12,:) = [62,1.92];
volt1_1(13,:) = [57,1.94];
volt1_1(14,:) = [53,2.00];
volt1_1(15,:) = [50,2.04];
volt1_1(16,:) = [47,2.13];

% volt1_1(17,:) = [44,2.09];
% volt1_1(18,:) = [42,2.15];
% volt1_1(19,:) = [40,2.20];



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
volt1_2(1,:) = [100000000000000,1.41];
volt1_2(2,:) = [272,1.50];
volt1_2(3,:) = [203,1.53];
volt1_2(4,:) = [162,1.57];
volt1_2(5,:) = [134,1.61];
volt1_2(6,:) = [115,1.65];

volt1_2(7,:) = [100,1.70];
volt1_2(8,:) = [89,1.74];
volt1_2(9,:) = [80,1.79];
volt1_2(10,:) = [73,1.84];
volt1_2(11,:) = [67,1.89];

volt1_2(12,:) = [62,1.92];
volt1_2(13,:) = [57,1.96];
volt1_2(14,:) = [53,1.99];
volt1_2(15,:) = [50,2.04];
volt1_2(16,:) = [47,2.06];

% volt1_2(17,:) = [44,2.09];
% volt1_2(18,:) = [42,2.14];
% volt1_2(19,:) = [40,2.19];



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
coff_sensor3_length=1/2*(coff_1+coff_2);
save('coff_sensor3_length.mat','coff_sensor3_length');