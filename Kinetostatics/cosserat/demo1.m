% �Ա�����ֽⷨ��Cosserat��

% 21-04-18
% ��һЩ����£����ַ�����õ���״����ͬ�ģ���ʹ��״��ͬ����õ�ĩ�˵����Ĵ�СҲ��ͬ�������Ի�
% ����һЩ����£����ַ�����õ���״��ͬ

% 21-04-20
% ������֮ǰ��bug���������ַ�����һ���ԽϺ�
% ������Ȼ�����������һ�µ���������� [0.6*L0;0.6*L0;0] ������״��ͬ����������λ�ò�ͬ

% ��Ȼ��һ��bugûŪ����
% poe tri cosserat ���ַ����������ĩ�˵����Ĵ�С����ͬ�ģ�����ĩ�����صĴ�С������ͬ
% ������������

% ����Ҫע����ǣ�dirac delta function���׸���ô����Ŀǰ�Ĵ���ʽ�ܼ�ª


% 21-04-20
% �޸���ĩ�����ص�bug ����tri��cosserat�����ĩ�˵��������ض���һ�µ�
% ��Ȼ��û����֤�����ǲ���poe�����ĩ�����ز�ͬ��ԭ���ǣ�
% poe����ʱʹ�õ��ǿռ��ſɱȣ�����ʱ��Ч�����õ���ȫ������ϵ��ԭ��
% ��tri��cosserat����ʱ��Ч�����õ���ĩ�ˣ�������߲�һ��ĩ�˾�ʸ�����

% 21-04-20��
% �ƺ�ĩ�˵�point force��point moment����Ҫ������f(s)��l(s)��
% �����Ѿ���ĩ�˵ı߽�������������
% ����������������Ժ�����ֽ���cosserat��õĽ������ȫһ���ˡ�

% cosserat��Ȼ����bug
% [0;0;pi]

% 21-04-21
% ֮ǰ��[0;0;pi]�������bug����ֵ����ѡ�õ��㷨�й�
% ���ʹ��Ĭ�ϵ�trust-region-dogleg�������zeros(6,1)Ϊ��ֵ�����޷���������
% ��Ҫ����ֵ��������ʵֵ���ӽ���ֵ���������
% ͬʱ�����ʹ��trust-region-dogleg������[0.3*L0;0.3*L0;pi/2]������������ֽⲻͬ�Ľ��

% ����[0;0;pi] ���ʹ��levenberg-marquardt��trust-region������zeros(6,1)Ϊ��ֵʱ��Ȼ�ܹ�������
% ����[0.3*L0;0.3*L0;pi/2] ���ʹ��levenberg-marquardt��trust-region�������������ֽ���ͬ�Ľ��

% Ŀǰ�Ľ����ǣ���Ҫ��Ĭ�ϵ�trust-region-dogleg

% ���⣬����ode45�Ĳ�����Ŀǰ�������趨����0.01 ʵ���ϲ���ȡ��ôС
% ���磬ȡ����Ϊ0.1�������Ч��Ҳ����

% 21-04-21
% ��������ֽⷨ����ʱҲ�����ʹ��trust-region-dogleg�㲻���������
% �����Ժ�ͱ���trust-region-dogleg�ˣ������ⲻ���ס�����

% 21-04-21
% ���԰�ı���������ڶ������
% ����[0.1*L0;0.3*L0;pi/3]
% cosseratʹ��3���㷨(LM TR TRD) �ֱ����������ģ̬�Ľ⣬���Ǻ�����˼��



%% ����ֽ�

clear
clc

wid=14e-3;
thi=0.3e-3;
E=197*1e9;

L0=1;
n=50;

pdes=[0;0;pi];
% pdes=[0.1*L0;0.3*L0;pi/3]
% pdes=[0.3*L0;0.3*L0;pi/2];
% pdes=[0.6*L0;0.6*L0;pi/2];

R1=planar_nR(E,L0,wid,thi,n,pdes);


f=@(x) cal_balance(x,R1);
options1 = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off','Algorithm','levenberg-marquardt');
x0=zeros(R1.n_seg+3,1);
[x_nR,fval_nR,exitflag_nR,output_nR] = fsolve(f,x0,options1);

R1.theta=x_nR(1:R1.n_seg);
R1.F=x_nR(R1.n_seg+1:R1.n_seg+3);

R1.update;
R1.plot_all;



%% Cosserat Model

wid=14e-3;
thi=0.3e-3;
E=197*1e9;

Iz=1/12*thi.^3*wid;

L=1;


f=@(x) check_balance(x,L,E,Iz,pdes);


x0=zeros(6,1);

options_a = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
options_b = optimoptions('fsolve','Display','off','Algorithm','trust-region');
options_c = optimoptions('fsolve','Display','off','Algorithm','trust-region-dogleg');
[x_cosserat,fval_cosserat,exitflag_cosserat,output_cosserat]=fsolve(f,x0,options_a);


% ��֤���

span = [0 L];
y0 = [0;0;0;x_cosserat(1:3)];
options3=odeset('MaxStep',1e-2);
[s,Y] = ode45(@(s,y) get_ydot(s,y,L,E,Iz,x_cosserat(4:6)), span, y0,options3);

hold on
plot(Y(:,1),Y(:,2),'r')

legend('principal axes decomposition','cosserat rod theory','Location','NorthWest')

%%
function res=check_balance(x,L,E,I,pdes)
    
    n0=x(1:2);
    m0=x(3);
    Fe=x(4:6);
    
    span = [0 L];
    y0 = [0;0;0;n0;m0];
    
    options=odeset('MaxStep',1e-2);
    
    [~,Y] = ode45(@(s,y) get_ydot(s,y,L,E,I,Fe), span, y0,options);
    
    ye=transpose(Y(end,:));
    
    res=zeros(6,1);
    res(1:2)=ye(1:2)-pdes(1:2);
    res(3)=ye(3)-pdes(3);
    res(4:5)=ye(4:5)-Fe(1:2);
    res(6)=transpose(ye(1:2))*[0 1;-1 0]*ye(4:5)+ye(6)-transpose(ye(1:2))*[0 1;-1 0]*Fe(1:2)-Fe(3);
    
end


function ydot=get_ydot(s,y,L,E,I,Fe)
    ydot=zeros(size(y));
    
    r=y(1:2);
    theta=y(3);
    n=y(4:5);
    m=y(6);
    
%     delta=1e-16;
%     if (L-s)<delta
%         f=Fe(1:2);
%         l=Fe(3);
%     else
%         f=[0;0];
%         l=0;
%     end
    
    f=[0;0];
    l=0;

    ydot(1:2)=[cos(theta);sin(theta)];
    ydot(3)=1/(E*I)*m;
    ydot(4:5)=-f;
    ydot(6)=-[cos(theta) sin(theta)]*[0 1;-1 0]*n-l;
end

function [r,J]=cal_balance(x,Rod)
    Rod.theta=x(1:Rod.n_seg);
    Rod.F=x(Rod.n_seg+1:Rod.n_seg+3);
    
    Rod.update;
    r=Rod.cal_r;
    J=Rod.cal_rdot;
end