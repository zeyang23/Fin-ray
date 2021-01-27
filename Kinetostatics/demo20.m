% ������ѧ��������

clear all;
clc;
addpath(genpath('.'));

%���ϵ��������
wid=5e-3;
thi=1e-3;
I=1/12*thi.^3*wid;
E=197*1e9;

%΢Ԫ����
n=50;

%�˵ĳ���
L0=0.2;

delta=L0/n;

%Ŀ��λ��
theta=60; % �Ƕ���

N=20;

end_pos=[0.6*L0;0.6*L0];

JACOB=@(x) Jacob_constraint_simple2(L0,I,E,n,x);

x_all=interp_solver3(theta,end_pos,N,n,L0,I,E,JACOB);

JACOB=@(x) Jacob_constraint_simple2(L0,I,E,n,x);


%��ʼ��̬����
g0=[eye(3),[L0;0;0];0,0,0,1];

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


figure
for i=1:N
    theta_solve=x_all(i,1:n);
    plot_pos2(w_all,q_all,g0,delta,theta_solve,L0)

end

rmpath(genpath('.'));