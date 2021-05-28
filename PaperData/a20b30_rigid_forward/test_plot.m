% 画图


%%
clear
clc

load('Nodes.mat');
load('Nodes_cal.mat');

Delta=3:3:18;

figure
hold on

subplot(4,2,1)
hold on
plot(Delta,Nodes.node1(:,1),'-o')
plot(Delta,Nodes_cal.node1_cal(:,1),'-*')
title("节点1",'Color','blue','FontSize',14 )
xlabel('夹钳水平移动距离/mm')
ylabel('节点x坐标/mm')
legend('实验测量','理论计算','Location','NorthWest')


subplot(4,2,2)
hold on
plot(Delta,Nodes.node5(:,1),'-o')
plot(Delta,Nodes_cal.node5_cal(:,1),'-*')
title("节点5",'Color','blue','FontSize',14)
xlabel('夹钳水平移动距离/mm')
ylabel('节点x坐标/mm')
legend('实验测量','理论计算','Location','NorthWest')


subplot(4,2,3)
hold on
plot(Delta,Nodes.node2(:,1),'-o')
plot(Delta,Nodes_cal.node2_cal(:,1),'-*')
title("节点2",'Color','blue','FontSize',14)
xlabel('夹钳水平移动距离/mm')
ylabel('节点x坐标/mm')
legend('实验测量','理论计算','Location','NorthWest')


subplot(4,2,4)
hold on
plot(Delta,Nodes.node6(:,1),'-o')
plot(Delta,Nodes_cal.node6_cal(:,1),'-*')
title("节点6",'Color','blue','FontSize',14)
xlabel('夹钳水平移动距离/mm')
ylabel('节点x坐标/mm')
legend('实验测量','理论计算','Location','NorthWest')


subplot(4,2,5)
hold on
plot(Delta,Nodes.node3(:,1),'-o')
plot(Delta,Nodes_cal.node3_cal(:,1),'-*')
title("节点3",'Color','blue','FontSize',14)
xlabel('夹钳水平移动距离/mm')
ylabel('节点x坐标/mm')
legend('实验测量','理论计算','Location','NorthWest')


subplot(4,2,6)
hold on
plot(Delta,Nodes.node7(:,1),'-o')
plot(Delta,Nodes_cal.node7_cal(:,1),'-*')
title("节点7",'Color','blue','FontSize',14)
xlabel('夹钳水平移动距离/mm')
ylabel('节点x坐标/mm')
legend('实验测量','理论计算','Location','NorthWest')

subplot(4,2,7)
hold on
plot(Delta,Nodes.node4(:,1),'-o')
plot(Delta,Nodes_cal.node4_cal(:,1),'-*')
title("节点4",'Color','blue','FontSize',14)
xlabel('夹钳水平移动距离/mm')
ylabel('节点x坐标/mm')
legend('实验测量','理论计算','Location','NorthWest')


subplot(4,2,8)
hold on
plot(Delta,Nodes.node8(:,1),'-o')
plot(Delta,Nodes_cal.node8_cal(:,1),'-*')
title("节点8",'Color','blue','FontSize',14)
xlabel('夹钳水平移动距离/mm')
ylabel('节点x坐标/mm')
legend('实验测量','理论计算','Location','NorthWest')


set(gcf,'Position',[100 100 800 800])

print('x_all','-djpeg','-r1500')




%%
clear
clc

load('Nodes.mat');
load('Nodes_cal.mat');

Delta=3:3:18;

figure
hold on

subplot(4,2,1)
hold on
plot(Delta,Nodes.node1(:,2),'-o')
plot(Delta,Nodes_cal.node1_cal(:,2),'-*')
title("节点1",'Color','blue','FontSize',14 )
xlabel('夹钳水平移动距离/mm')
ylabel('节点y坐标/mm')
legend('实验测量','理论计算','Location','NorthWest')

subplot(4,2,2)
hold on
plot(Delta,Nodes.node5(:,2),'-o')
plot(Delta,Nodes_cal.node5_cal(:,2),'-*')
title("节点5",'Color','blue','FontSize',14)
xlabel('夹钳水平移动距离/mm')
ylabel('节点y坐标/mm')
legend('实验测量','理论计算','Location','NorthWest')

subplot(4,2,3)
hold on
plot(Delta,Nodes.node2(:,2),'-o')
plot(Delta,Nodes_cal.node2_cal(:,2),'-*')
title("节点2",'Color','blue','FontSize',14)
xlabel('夹钳水平移动距离/mm')
ylabel('节点y坐标/mm')
legend('实验测量','理论计算','Location','NorthWest')

subplot(4,2,4)
hold on
plot(Delta,Nodes.node6(:,2),'-o')
plot(Delta,Nodes_cal.node6_cal(:,2),'-*')
title("节点6",'Color','blue','FontSize',14)
xlabel('夹钳水平移动距离/mm')
ylabel('节点y坐标/mm')
legend('实验测量','理论计算','Location','NorthWest')

subplot(4,2,5)
hold on
plot(Delta,Nodes.node3(:,2),'-o')
plot(Delta,Nodes_cal.node3_cal(:,2),'-*')
title("节点3",'Color','blue','FontSize',14)
xlabel('夹钳水平移动距离/mm')
ylabel('节点y坐标/mm')
legend('实验测量','理论计算','Location','NorthWest')

subplot(4,2,6)
hold on
plot(Delta,Nodes.node7(:,2),'-o')
plot(Delta,Nodes_cal.node7_cal(:,2),'-*')
title("节点7",'Color','blue','FontSize',14)
xlabel('夹钳水平移动距离/mm')
ylabel('节点y坐标/mm')
legend('实验测量','理论计算','Location','NorthWest')

subplot(4,2,7)
hold on
plot(Delta,Nodes.node4(:,2),'-o')
plot(Delta,Nodes_cal.node4_cal(:,2),'-*')
title("节点4",'Color','blue','FontSize',14)
xlabel('夹钳水平移动距离/mm')
ylabel('节点y坐标/mm')
legend('实验测量','理论计算','Location','NorthWest')

subplot(4,2,8)
hold on
plot(Delta,Nodes.node8(:,2),'-o')
plot(Delta,Nodes_cal.node8_cal(:,2),'-*')
title("节点8",'Color','blue','FontSize',14)
xlabel('夹钳水平移动距离/mm')
ylabel('节点y坐标/mm')
legend('实验测量','理论计算','Location','NorthWest')


set(gcf,'Position',[100 100 800 800])

print('y_all','-djpeg','-r1500')