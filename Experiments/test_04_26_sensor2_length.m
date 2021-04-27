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
volt1_1(1,:) = [100000000000000,3.09];
volt1_1(2,:) = [272,3.18];
volt1_1(3,:) = [203,3.24];
volt1_1(4,:) = [162,3.31];
volt1_1(5,:) = [134,3.37];
volt1_1(6,:) = [115,3.42];

volt1_1(7,:) = [100,3.46];
volt1_1(8,:) = [89,3.52];
volt1_1(9,:) = [80,3.60];
volt1_1(10,:) = [73,3.68];
volt1_1(11,:) = [67,3.74];

volt1_1(12,:) = [62,3.80];
volt1_1(13,:) = [57,3.87];
volt1_1(14,:) = [53,3.94];
volt1_1(15,:) = [50,4.01];
volt1_1(16,:) = [47,4.06];

% volt1_1(17,:) = [44,4.12];
% volt1_1(18,:) = [42,4.21];
% volt1_1(19,:) = [40,4.30];


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
volt1_2(1,:) = [100000000000000,3.09];
volt1_2(2,:) = [272,3.17];
volt1_2(3,:) = [203,3.22];
volt1_2(4,:) = [162,3.28];
volt1_2(5,:) = [134,3.35];
volt1_2(6,:) = [115,3.40];

volt1_2(7,:) = [100,3.48];
volt1_2(8,:) = [89,3.54];
volt1_2(9,:) = [80,3.61];
volt1_2(10,:) = [73,3.67];
volt1_2(11,:) = [67,3.75];

volt1_2(12,:) = [62,3.81];
volt1_2(13,:) = [57,3.86];
volt1_2(14,:) = [53,3.93];
volt1_2(15,:) = [50,3.97];
volt1_2(16,:) = [47,4.04];

% volt1_2(17,:) = [44,4.14];
% volt1_2(18,:) = [42,4.25];
% volt1_2(19,:) = [40,4.32];



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
coff_sensor2_length=1/2*(coff_1+coff_2);
save('coff_sensor2_length.mat','coff_sensor2_length');