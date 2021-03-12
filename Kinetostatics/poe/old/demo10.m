% ��֤���������ĸ�
% ʹ��ţ�ٷ� �����ſɱȾ���
% FȡΪ��ֵ������������
% n��΢Ԫ������ȡn����������ͬ��demo8
% F�ķ������⣬һֱ�����

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
K=diag((1/c)*ones(n,1));

%����
F=[0;-20;0;0;0;0];

%��ʼ��̬����
g0=[eye(3),[L;0;0];0,0,0,1];

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


%���ƽ��λ��
f=@(theta) cal_balance(K,w_all,q_all,g0,theta,F);
JACOB=@(theta) Jacob_solve(K,w_all,q_all,theta,F);

theta0=linspace(0,pi/100,n);
TOL=1e-6;

theta_solve=Newton_nd(f,JACOB,theta0,TOL);

theta_fsolve=fsolve(f,linspace(0,pi/100,n));

error=norm(theta_solve-theta_fsolve);

clc
disp("difference between Newton and fsolve")
disp(error)

[g_exp,~]=exp_fkine(w_all,q_all,g0,theta_solve);

%�õ������˵�λ��
plot_pos(w_all,q_all,g0,delta,theta_solve)

save('G1.mat','g_exp');
save('Theta1.mat','theta_solve');
save('F1.mat','F');

rmpath(genpath('.'));

function J=Jacob_solve(K,w_all,q_all,theta_all,F)
    L=0.1;
    g0=[eye(3),[L;0;0];0,0,0,1];
    [Jacobs,~,~]=exp_jacob(w_all,q_all,g0,theta_all);
    J=K-partial_J_theta(w_all,q_all,theta_all,Jacobs,F);
end