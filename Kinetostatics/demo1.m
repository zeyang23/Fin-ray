% ����ָ������������˶�ѧ

clear all;
clc;
addpath(genpath('.'));


%��ʼ��̬����
g0=[eye(3),[2;0;2];0,0,0,1];  

%���е�w�������ɾ���
w_all=[[0;0;1],[0;1;0],[1/sqrt(5);0;2/sqrt(5)]]; 

%���е�q�������ɾ���
q_all=[[0;0;0],[0;0;2],[0;0;0]];

theta_all=[pi/10;-pi/6;pi/12];   

%ָ���������㣬����w_all��q_all��g0��theta_all
[g_st,xi_all]=exp_fkine(w_all,q_all,g0,theta_all);  

rmpath(genpath('.'));