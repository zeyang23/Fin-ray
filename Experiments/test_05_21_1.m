% 标定弯曲传感器

clear
clc

%% 标定第I块弯曲传感器
% 测试条件：电流4mA
volt1 = [];
volt1(1,:) = [100000000000000,2.43];
volt1(2,:) = [272,2.53];
volt1(3,:) = [203,2.57];
volt1(4,:) = [162,2.62];
volt1(5,:) = [134,2.67];
volt1(6,:) = [115,2.72];

volt1(7,:) = [100,2.78];
volt1(8,:) = [89,2.82];
volt1(9,:) = [80,2.89];
volt1(10,:) = [73,2.93];
volt1(11,:) = [67,2.99];


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
volt2(1,:) = [100000000000000,2.17];
volt2(2,:) = [272,2.27];
volt2(3,:) = [203,2.31];
volt2(4,:) = [162,2.345];
volt2(5,:) = [134,2.38];
volt2(6,:) = [115,2.42];

volt2(7,:) = [100,2.46];
volt2(8,:) = [89,2.50];
volt2(9,:) = [80,2.55];
volt2(10,:) = [73,2.60];
volt2(11,:) = [67,2.65];


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



%% 标定第III块弯曲传感器
volt3 = [];
volt3(1,:) = [100000000000000,2.06];
volt3(2,:) = [272,2.14];
volt3(3,:) = [203,2.17];
volt3(4,:) = [162,2.21];
volt3(5,:) = [134,2.25];
volt3(6,:) = [115,2.30];

volt3(7,:) = [100,2.34];
volt3(8,:) = [89,2.38];
volt3(9,:) = [80,2.42];
volt3(10,:) = [73,2.46];
volt3(11,:) = [67,2.50];


% 减掉初值
init3=volt3(1,2);
volt3(:,2)=volt3(:,2)-init3;


volt3(:,1) = 1./volt3(:,1)*38;
plot(volt3(:,2),volt3(:,1),'ob','linewidth',1);

coff_3 = polyfit(volt3(:,2),volt3(:,1),1);

hold on
grid on
y_3 = polyval(coff_3,volt3(:,2));
plot(volt3(:,2),y_3,'--b','linewidth',1)