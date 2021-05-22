% 标定弯曲传感器

clear
clc

%% 标定第I块弯曲传感器
% 测试条件：电流5.6mA
volt1 = [];
volt1(1,:) = [100000000000000,1.58];
volt1(2,:) = [272,1.62];
volt1(3,:) = [203,1.68];
volt1(4,:) = [162,1.72];
volt1(5,:) = [134,1.76];
volt1(6,:) = [115,1.81];

volt1(7,:) = [100,1.85];
volt1(8,:) = [89,1.89];
volt1(9,:) = [80,1.94];
volt1(10,:) = [73,1.97];
volt1(11,:) = [67,2.02];


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


%% 标定第II块弯曲传感器
volt2 = [];
volt2(1,:) = [100000000000000,1.52];
volt2(2,:) = [272,1.61];
volt2(3,:) = [203,1.65];
volt2(4,:) = [162,1.69];
volt2(5,:) = [134,1.73];
volt2(6,:) = [115,1.78];

volt2(7,:) = [100,1.83];
volt2(8,:) = [89,1.87];
volt2(9,:) = [80,1.91];
volt2(10,:) = [73,1.95];
volt2(11,:) = [67,1.995];


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

