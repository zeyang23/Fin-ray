%��֤�˵����

clear all;
clc;
addpath(genpath('.'));

%��ɢ��΢Ԫ����
n=50;

%���Ը˵��ܳ���
L=1;

%΢Ԫ�ĳ���
delta=L/n;

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

theta_all=pi/100*ones(1,n);
theta_all=[theta_all,0];

[g_exp,~]=exp_fkine(w_all,q_all,g0,theta_all);


%�õ������˵�λ��
pos=[0;0];
for i =1:n
    g0i=[eye(3),[(i-0.5)*delta;0;0];0,0,0,1];
    [gsti,~]=exp_fkine(w_all(:,1:i),q_all(:,1:i),g0i,theta_all(:,1:i));
    pos=[pos,gsti(1:2,4)];
end

pos=[pos,g_exp(1:2,4)];

plot(pos(1,:),pos(2,:),'-o','MarkerSize',2)

rmpath(genpath('.'));