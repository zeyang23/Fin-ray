% Լ���µĸ� �̶��˳�
% ʹ��ţ�ٷ����������ݶ�
% �𲽱ƽ�Ŀ��λ��

clear all;
clc;
addpath(genpath('.'));

%���ϵ��������
I=0.25*pi*1e-12;
E=56*1e9;

%΢Ԫ����
n=10;

%�˵ĳ���
L0=1;

delta=L0/n;

%����ĩ��λ��
theta=45;
end_pos=[0.8*L0;0.6*L0];

JACOB=@(x) Jacob_constraint_simple(L0,I,E,n,x);

xsolve=interp_solver(theta,end_pos,10,n,L0,I,E,JACOB);

theta_solve=(xsolve(1:end-6))';


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

[g_exp,~]=exp_fkine(w_all,q_all,g0,theta_solve);

%�õ������˵�λ��
plot_pos(w_all,q_all,g0,delta,theta_solve)

rmpath(genpath('.'));