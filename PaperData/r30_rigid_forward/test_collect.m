%% 汇总6组数据
clear
clc

Nodes=struct;

Nodes.node1=[];
Nodes.node2=[];
Nodes.node3=[];
Nodes.node4=[];
Nodes.node5=[];
Nodes.node6=[];
Nodes.node7=[];
Nodes.node8=[];


Nodes_cal=struct;

Nodes_cal.node1_cal=[];
Nodes_cal.node2_cal=[];
Nodes_cal.node3_cal=[];
Nodes_cal.node4_cal=[];
Nodes_cal.node5_cal=[];
Nodes_cal.node6_cal=[];
Nodes_cal.node7_cal=[];
Nodes_cal.node8_cal=[];

%%

i=6;

load('r30_delta30.mat')
load('r30_delta30_cal.mat')

Nodes.node1(i,:)=loc_points_real(1,:);
Nodes.node2(i,:)=loc_points_real(2,:);
Nodes.node3(i,:)=loc_points_real(3,:);
Nodes.node4(i,:)=loc_points_real(4,:);
Nodes.node5(i,:)=loc_points_real(5,:);
Nodes.node6(i,:)=loc_points_real(6,:);
Nodes.node7(i,:)=loc_points_real(7,:);
Nodes.node8(i,:)=loc_points_real(8,:);

Nodes_cal.node1_cal(i,:)=rigid_abs_pos(1,:);
Nodes_cal.node2_cal(i,:)=rigid_abs_pos(2,:);
Nodes_cal.node3_cal(i,:)=rigid_abs_pos(3,:);
Nodes_cal.node4_cal(i,:)=rigid_abs_pos(4,:);
Nodes_cal.node5_cal(i,:)=rigid_abs_pos(5,:);
Nodes_cal.node6_cal(i,:)=rigid_abs_pos(6,:);
Nodes_cal.node7_cal(i,:)=rigid_abs_pos(7,:);
Nodes_cal.node8_cal(i,:)=rigid_abs_pos(8,:);