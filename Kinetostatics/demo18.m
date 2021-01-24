% Լ���µĸ� �̶��˳�
% ʹ��ţ�ٷ����������ݶ�
% ʹ�������ſɱ� Jb

% rotz����������ǽǶȣ����ǻ��ȣ�
% ����Ч���ͺܺ���

% ��������λ�ˣ�ʵ�ֶ���
% ����״̬��ʼ��������λ�˵ķ��������⡣Ŀǰ�ǵ����㡣

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
theta=45; % �Ƕ���

N=50;

end_pos=[0.5*L0;0.5*L0];

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

avi_object = VideoWriter('demo18.avi');
avi_object.FrameRate = 10;
open(avi_object);
figure
for I=1:N
    theta_solve=x_all(I,1:n);
    plot_pos2(w_all,q_all,g0,delta,theta_solve,L0)

    M = getframe;
    writeVideo(avi_object,M);

    if I < N
        clf;
    end
end
close(avi_object);

rmpath(genpath('.'));