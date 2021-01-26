% Լ���µĸ� �̶��˳�
% ʹ��ţ�ٷ����������ݶ�
% ʹ�������ſɱ� Jb

% rotz����������ǽǶȣ����ǻ��ȣ�
% ����Ч���ͺܺ���

% ����ʹ�ÿռ��ſɱ�

% �������⣺���һ���ؽڵ�����Ǵ��

clear all;
clc;
addpath(genpath('.'));

%���ϵ��������
wid=5e-3;
thi=1e-3;
I=1/12*thi.^3*wid;
E=197*1e9;

%΢Ԫ����
n=20;

%�˵ĳ���
L0=0.2;

delta=L0/n;

%Ŀ��λ��
theta=60; % �Ƕ���

%Ŀ��λ��
gt=[rotz(theta),[0.6*L0;0.6*L0;0];0,0,0,1];

f=@(x) cal_constraint_simple(L0,I,E,n,gt,x);

x0=zeros(n+6,1);

x0(1:end-6)=linspace(0,0,n);

x0(end-5:end)=0*ones(6,1);

% load('x_demo13.mat')
% x0=xsolve;

JACOB=@(x) Jacob_constraint_simple(L0,I,E,n,x);

TOL=1e-3;

% x_fsolve=fsolve(f,x0);

xsolve=Newton_nd(f,JACOB,x0,TOL);

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