%用RTB校验指数坐标算法的正确性
clear all;
clc;
addpath(genpath('.'));

%% 指数坐标法
%离散的微元个数
n=8;

%柔性杆的总长度
L=1;

%微元的长度
delta=L/n;

%初始姿态矩阵
g0=[eye(3),[L;0;0];0,0,0,1];

%所有的w向量构成矩阵
w_all=[];
for i=1:n
    w_all=[w_all,[0;0;1]];
end
w_all=[w_all,[0;0;1]];

%所有的q向量构成矩阵
q_all=zeros(3,n);
for i=1:n
    q_all(:,i)=[1;0;0]*delta*(i-1/2);
end
q_all=[q_all,[1;0;0]*L];

theta_all=pi/100*ones(1,n);
theta_all=[theta_all,0];

[g_exp,~]=exp_fkine(w_all,q_all,g0,theta_all);


%得到整个杆的位姿
pos=[0;0];
for i =1:n
    g0i=[eye(3),[(i-0.5)*delta;0;0];0,0,0,1];
    [gsti,~]=exp_fkine(w_all(:,1:i),q_all(:,1:i),g0i,theta_all(:,1:i));
    pos=[pos,gsti(1:2,4)];
end

pos=[pos,g_exp(1:2,4)];

plot(pos(1,:),pos(2,:),'-*')

%% RTB
% D-H table
a1 = 0.5*delta; alpha1 = 0; d1 = 0; theta1 = 0;
a2 = delta; alpha2 = 0; d2 = 0; theta2 = 0;
a3 = delta; alpha3 = 0; d3 = 0; theta3 = 0;
a4 = delta; alpha4 = 0; d4 = 0; theta4 = 0;
a5 = delta; alpha5 = 0; d5 = 0; theta5 = 0;
a6 = delta; alpha6 = 0; d6 = 0; theta6 = 0;
a7 = delta; alpha7 = 0; d7 = 0; theta7 = 0;
a8 = delta; alpha8 = 0; d8 = 0; theta8 = 0;
a9 = 0.5*delta; alpha9 = 0; d9 = 0; theta9 = 0;
% a10 =delta; alpha10 = 0;d10 = 0; theta10 = 0;
% a11 =0.5*delta;alpha11 = 0;d11 = 0; theta11 = 0;

%generate links
LINK(1) = Link('d',d1,'a',a1,'alpha',alpha1);
LINK(2) = Link('d',d2,'a',a2,'alpha',alpha2);
LINK(3) = Link('d',d3,'a',a3,'alpha',alpha3);
LINK(4) = Link('d',d4,'a',a4,'alpha',alpha4);
LINK(5) = Link('d',d5,'a',a5,'alpha',alpha5);
LINK(6) = Link('d',d6,'a',a6,'alpha',alpha6);
LINK(7) = Link('d',d7,'a',a7,'alpha',alpha7);
LINK(8) = Link('d',d8,'a',a8,'alpha',alpha8);
LINK(9) = Link('d',d9,'a',a9,'alpha',alpha9);

%generate robotic
rod = SerialLink(LINK, 'name', 'rod');

Q=[0,theta_all];
Q(end)=[];

%plot
figure
rod.plot(Q,'view',[0,90]);
g_rtb=rod.fkine(Q);

jacobs_rtb=rod.jacob0(Q);
jacobe_rtb=rod.jacobe(Q);

[jacobs_exp,jacobe_exp,~]=exp_jacob(w_all,q_all,g0,theta_all);


g_rtb
g_exp

jacobe_rtb
jacobe_exp

rmpath(genpath('.'));