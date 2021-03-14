% Լ���µĸˣ�ʹ���������ݶ��Ż��㷨

clear all;
clc;
addpath(genpath('.'));

%���ϵ���������
I=0.25*pi*1e-12;
E=56*1e9;

%΢Ԫ����
n=5;

%��ʼ����
L0=0.1;

%Ŀ��λ��
gt=[eye(3),[0.05;0.05;0];0,0,0,1];

f=@(x) cal_constraint(I,E,n,gt,x);

x0=zeros(n+7,1);
x0(1)=L0;
x0(8:end)=linspace(0,0,n);

xsolve=fsolve(f,x0);

theta_final=(xsolve(2:end-6))';
L_final=xsolve(1);
delta=L_final/n;

%�õ������˵�λ��

%���е�w�������ɾ���
w_all=[];
for i=1:n
    w_all=[w_all,[0;0;1]];
end

%���е�q�������ɾ���
q_all=zeros(3,n);
for i=1:n
    q_all(:,i)=[1;0;0]*delta*(i-1/2);
end

pos=[0;0];
g0=[eye(3),[L_final;0;0];0,0,0,1];
[g_exp,~]=exp_fkine(w_all,q_all,g0,theta_final);
for i =1:n
    g0i=[eye(3),[(i-0.5)*delta;0;0];0,0,0,1];
    [gsti,~]=exp_fkine(w_all(:,1:i-1),q_all(:,1:i-1),g0i,theta_final(:,1:i-1));
    pos=[pos,gsti(1:2,4)];
end
pos=[pos,g_exp(1:2,4)];

plot(pos(1,:),pos(2,:),'-o')

rmpath(genpath('.'));