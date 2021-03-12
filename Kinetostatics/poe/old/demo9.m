% Լ���µĸ�
% ʹ��ţ�ٷ����������ݶ�

clear all;
clc;
addpath(genpath('.'));

%���ϵ��������
I=0.25*pi*1e-12;
E=56*1e9;

%΢Ԫ����
n=6;

%��ʼ����
L0=0.1;

%Ŀ��λ��
gt=[eye(3),[0.1;0;0];0,0,0,1];

f=@(x) cal_constraint(I,E,n,gt,x);


x0=zeros(n+7,1);
x0(1)=L0;
x0(2:end-6)=linspace(0,pi/10,n);

JACOB=@(x) Jacob_constraint(I,E,n,x);

TOL=1e-6;

xsolve=Newton_nd(f,JACOB,x0,TOL);



theta_final=(xsolve(8:end))';
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

plot_pos(w_all,q_all,g0,delta,theta_solve)

rmpath(genpath('.'));