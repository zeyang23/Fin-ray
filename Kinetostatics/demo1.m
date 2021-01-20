% 基于指数坐标的正向运动学

clear all;
clc;
addpath(genpath('.'));


%初始姿态矩阵
g0=[eye(3),[2;0;2];0,0,0,1];  

%所有的w向量构成矩阵
w_all=[[0;0;1],[0;1;0],[1/sqrt(5);0;2/sqrt(5)]]; 

%所有的q向量构成矩阵
q_all=[[0;0;0],[0;0;2],[0;0;0]];

theta_all=[pi/10;-pi/6;pi/12];   

%指数正向运算，输入w_all，q_all，g0，theta_all
[g_st,xi_all]=exp_fkine(w_all,q_all,g0,theta_all);  

rmpath(genpath('.'));