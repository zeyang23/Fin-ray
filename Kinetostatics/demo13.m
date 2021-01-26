% Լ���µĸ� �̶��˳�
% ʹ��ţ�ٷ����������ݶ�
% ʹ�������ſɱ� Jb

% rotz����������ǽǶȣ����ǻ��ȣ�
% ����Ч���ͺܺ���

clear all;
clc;
addpath(genpath('.'));

%���ϵ��������
I=0.25*pi*1e-12;
E=56*1e9;

%΢Ԫ����
n=50;

%�˵ĳ���
L0=1;

delta=L0/n;

%Ŀ��λ��
gt=[rotz(100),[0.5*L0;0.5*L0;0];0,0,0,1];

f=@(x) cal_constraint_simple2(L0,I,E,n,gt,x);

x0=zeros(n+6,1);

x0(1:end-6)=linspace(0,0,n);

x0(end-5:end)=0*ones(6,1);

JACOB=@(x) Jacob_constraint_simple2(L0,I,E,n,x);

TOL=1e-3;

% x_fsolve=fsolve(f,x0);

xsolve=Newton_nd(f,JACOB,x0,TOL);
save('x_demo13.mat','xsolve')

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