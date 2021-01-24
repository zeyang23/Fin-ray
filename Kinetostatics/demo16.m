% Լ���µĸ� �̶��˳�
% ʹ��ţ�ٷ����������ݶ�
% ʹ�������µĸ˵Ľ����У���ſɱȵ���ȷ�� ����demo10
% ���ʹ�������ſɱȣ�Ч���ȿռ��ſɱȺúܶ�
% ���⣺����Լ����ָ��ĩ����̬��λ�ã������㲻������������Ǵ������ĸ����ɵ���ʵ���
% �����𲽱ƽ�Ŀ��λ��

clear all;
clc;
addpath(genpath('.'));

%���ϵ��������
I=0.25*pi*1e-12;
E=56*1e9;

%΢Ԫ����
n=50;

%�˵ĳ���
L0=0.1;

delta=L0/n;

theta=180;
end_pos=[0.5*L0;0.5*L0];

JACOB=@(x) Jacob_constraint_simple2(L0,I,E,n,x);

xsolve=interp_solver2(theta,end_pos,10,n,L0,I,E,JACOB);

TOL=1e-3;

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

[g_solve,~]=exp_fkine(w_all,q_all,g0,theta_solve);
disp('g_solve')
disp(g_solve)

%�õ������˵�λ��
plot_pos(w_all,q_all,g0,delta,theta_solve)

rmpath(genpath('.'));