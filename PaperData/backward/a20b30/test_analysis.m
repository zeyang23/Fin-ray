clear
clc

load('Nodes.mat');
load('Nodes_cal.mat');

ex1=sum(abs(Nodes.node1(:,1)-Nodes_cal.node1_cal(:,1)))/6;
ex2=sum(abs(Nodes.node2(:,1)-Nodes_cal.node2_cal(:,1)))/6;
ex3=sum(abs(Nodes.node3(:,1)-Nodes_cal.node3_cal(:,1)))/6;
ex4=sum(abs(Nodes.node4(:,1)-Nodes_cal.node4_cal(:,1)))/6;
ex5=sum(abs(Nodes.node5(:,1)-Nodes_cal.node5_cal(:,1)))/6;
ex6=sum(abs(Nodes.node6(:,1)-Nodes_cal.node6_cal(:,1)))/6;
ex7=sum(abs(Nodes.node7(:,1)-Nodes_cal.node7_cal(:,1)))/6;
ex8=sum(abs(Nodes.node8(:,1)-Nodes_cal.node8_cal(:,1)))/6;

ex=1/8*(ex1+ex2+ex3+ex4+ex5+ex6+ex7+ex8);

ey1=sum(abs(Nodes.node1(:,2)-Nodes_cal.node1_cal(:,2)))/6;
ey2=sum(abs(Nodes.node2(:,2)-Nodes_cal.node2_cal(:,2)))/6;
ey3=sum(abs(Nodes.node3(:,2)-Nodes_cal.node3_cal(:,2)))/6;
ey4=sum(abs(Nodes.node4(:,2)-Nodes_cal.node4_cal(:,2)))/6;
ey5=sum(abs(Nodes.node5(:,2)-Nodes_cal.node5_cal(:,2)))/6;
ey6=sum(abs(Nodes.node6(:,2)-Nodes_cal.node6_cal(:,2)))/6;
ey7=sum(abs(Nodes.node7(:,2)-Nodes_cal.node7_cal(:,2)))/6;
ey8=sum(abs(Nodes.node8(:,2)-Nodes_cal.node8_cal(:,2)))/6;

ey=1/8*(ey1+ey2+ey3+ey4+ey5+ey6+ey7+ey8);
