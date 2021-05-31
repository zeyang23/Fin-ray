% 标定弯曲传感器

% 21-05-31
% 应该使用正比例函数拟合，不应该用一次函数拟合

clear
clc

%% 标定第I块弯曲传感器
% 测试条件：电流4mA
volt1 = [];
volt1(1,:) = [100000000000000,1.81];
volt1(2,:) = [272,1.94];
volt1(3,:) = [203,1.98];
volt1(4,:) = [162,2.03];
volt1(5,:) = [134,2.08];
volt1(6,:) = [115,2.13];

volt1(7,:) = [100,2.19];
volt1(8,:) = [89,2.24];
volt1(9,:) = [80,2.30];
volt1(10,:) = [73,2.34];
volt1(11,:) = [67,2.40];


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
volt2(1,:) = [100000000000000,2.09];
volt2(2,:) = [272,2.22];
volt2(3,:) = [203,2.275];
volt2(4,:) = [162,2.33];
volt2(5,:) = [134,2.37];
volt2(6,:) = [115,2.43];

volt2(7,:) = [100,2.50];
volt2(8,:) = [89,2.55];
volt2(9,:) = [80,2.63];
volt2(10,:) = [73,2.675];
volt2(11,:) = [67,2.73];


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



%% 标定第III块弯曲传感器
volt3 = [];
volt3(1,:) = [100000000000000,1.84];
volt3(2,:) = [272,1.96];
volt3(3,:) = [203,2.00];
volt3(4,:) = [162,2.06];
volt3(5,:) = [134,2.11];
volt3(6,:) = [115,2.16];

volt3(7,:) = [100,2.20];
volt3(8,:) = [89,2.25];
volt3(9,:) = [80,2.30];
volt3(10,:) = [73,2.35];
volt3(11,:) = [67,2.395];


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

coff_3_real=volt3(:,2)\volt3(:,1);