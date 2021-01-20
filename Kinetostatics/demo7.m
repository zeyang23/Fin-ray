% ��֤���������ĸ�
% ��ת�ſɱȾ���demo4�е�Jδ����ת�������Ǵ��
% demo7��demo4�õ��Ľ����ͬ
% y���������û������

clear all;
clc;
addpath(genpath('.'));

%��ɢ��΢Ԫ����
n=50;

%���Ը˵��ܳ���
L=0.1;

%΢Ԫ�ĳ���
delta=L/n;

%���ϵ��������
I=0.25*pi*1e-12;
E=56*1e9;
c=delta/(E*I);

%�նȾ���
K=diag((1/c)*ones(n+1,1));

%����
F=[0;0;0;0;-10;0];

%��ʼ��̬����
g0=[eye(3),[L;0;0];0,0,0,1];

%���е�w�������ɾ���
w_all=[];
for i=1:n
    w_all=[w_all,[0;0;1]];
end
w_all=[w_all,[0;0;1]];

%���е�q�������ɾ���
q_all=zeros(3,n);
for i=1:n
    q_all(:,i)=[1;0;0]*delta*(i-1/2);
end
q_all=[q_all,[1;0;0]*L];

%���ƽ��λ��
f=@(theta) cal_balance2(K,w_all,q_all,g0,theta,F);
theta_solve=fsolve(f,linspace(0,pi/100,n+1));

[g_exp,~]=exp_fkine(w_all,q_all,g0,theta_solve);

%�õ������˵�λ��
pos=[0;0];
for i =1:n
    g0i=[eye(3),[(i-0.5)*delta;0;0];0,0,0,1];
    [gsti,~]=exp_fkine(w_all(:,1:i),q_all(:,1:i),g0i,theta_solve(:,1:i));
    pos=[pos,gsti(1:2,4)];
end
pos=[pos,g_exp(1:2,4)];

plot(pos(1,:),pos(2,:),'-o','MarkerSize',2)

rmpath(genpath('.'));