% ����ָ����������ſɱȾ���

clear all;
clc;
addpath(genpath('.'));

%��ʼ��̬����
g0=[[0 0 1;0 -1 0;1 0 0],[1.6;0;0.5];0,0,0,1];

%���е�w�������ɾ���
w_all=[[0;0;1],[0;1;0],[0;1;0],[1;0;0],[0;1;0],[1;0;0]];

%���е�q�������ɾ���
q_all=[[0;0;0],[0;0;0.5],[0.5;0;0.5],[0;0;0.5],[1.1;0;0.5],[0;0;0.5]];

theta_all=pi/20*[1;2;3;4;5;6];   

[g_st,xi_all]=exp_fkine(w_all,q_all,g0,theta_all);   

[jacobs,jacobe,gap]=exp_jacob(w_all,q_all,g0,theta_all);

rmpath(genpath('.'));