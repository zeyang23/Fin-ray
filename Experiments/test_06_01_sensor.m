% 标定弯曲传感器

% 21-05-31
% 应该使用正比例函数拟合，不应该用一次函数拟合

clear
clc

%% 标定第I块弯曲传感器
% 测试条件：电流4mA
volt1 = [];
volt1(1,:) = [100000000000000,2.19];
volt1(2,:) = [272,2.33];
volt1(3,:) = [203,2.38];
volt1(4,:) = [162,2.445];
volt1(5,:) = [134,2.51];
volt1(6,:) = [115,2.57];

volt1(7,:) = [100,2.64];
volt1(8,:) = [89,2.685];
volt1(9,:) = [80,2.75];
volt1(10,:) = [73,2.80];
volt1(11,:) = [67,2.86];


% 减掉初值
init1=volt1(1,2);
volt1(:,2)=volt1(:,2)-init1;


volt1(:,1) = 1./volt1(:,1)*38;
plot(volt1(:,2),volt1(:,1),'or','linewidth',1);

coff_1 = polyfit(volt1(:,2),volt1(:,1),1);

hold on
grid on
y_1 = polyval(coff_1,volt1(:,2));
plot(volt1(:,2),y_1,'--r','linewidth',1)

coff_1_real=volt1(:,2)\volt1(:,1);


%% 标定第II块弯曲传感器
volt2 = [];
volt2(1,:) = [100000000000000,2.11];
volt2(2,:) = [272,2.27];
volt2(3,:) = [203,2.32];
volt2(4,:) = [162,2.38];
volt2(5,:) = [134,2.44];
volt2(6,:) = [115,2.50];

volt2(7,:) = [100,2.56];
volt2(8,:) = [89,2.61];
volt2(9,:) = [80,2.69];
volt2(10,:) = [73,2.73];
volt2(11,:) = [67,2.79];


% 减掉初值
init2=volt2(1,2);
volt2(:,2)=volt2(:,2)-init2;


volt2(:,1) = 1./volt2(:,1)*38;
plot(volt2(:,2),volt2(:,1),'og','linewidth',1);

coff_2 = polyfit(volt2(:,2),volt2(:,1),1);

hold on
grid on
y_2 = polyval(coff_2,volt2(:,2));
plot(volt2(:,2),y_2,'--g','linewidth',1)


coff_2_real=volt2(:,2)\volt2(:,1);