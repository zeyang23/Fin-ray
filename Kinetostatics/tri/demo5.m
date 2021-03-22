% ���ԣ�ͨ��3�������������ݷ������԰����״
% ʹ��ţ�ٷ����

% 21-03-22���Խ����Ч�������ף�����������������
% �������ɹ��������������״�����λ����ƫ�ƣ�������������״�����ף������󲻳�������������

% ��֧��ʹ��zeros(n+3,1)��Ϊ��ֵ��
% partial��Ӧ�ô�


% ����3��Delta
clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=40;

pdes=[0.6*L0;0.6*L0;2/3*pi];   % ���������Ч���ǳ���
% pdes=[0.6*L0;0.6*L0;1/3*pi];   % ��������������״ƫ���е�Զ
% pdes=[0.6*L0;0.6*L0;0];        % �������������
% pdes=[0.3*L0;0.3*L0;pi/2];     % �������Ҳ������

R1=planar_nR(E,L0,wid,thi,n,pdes);

TOL=1e-6;
R1.Newton(TOL);
R1.plot_all;
hold on

% �������Ĵ�С��λ��
sensor_length=1/8;

p1=fix((1/4-sensor_length/2)*n);
q1=fix((1/4+sensor_length/2)*n);

p2=fix((2/4-sensor_length/2)*n);
q2=fix((2/4+sensor_length/2)*n);

p3=fix((3/4-sensor_length/2)*n);
q3=fix((3/4+sensor_length/2)*n);

Delta1=R1.theta(q1)-R1.theta(p1);
Delta2=R1.theta(q2)-R1.theta(p2);
Delta3=R1.theta(q3)-R1.theta(p3);


% ���Է���
Rod=planar_nR(E,L0,wid,thi,n,[0;0;0]);

D1=zeros(1,n);
D1(p1)=-1;
D1(q1)=1;

D2=zeros(1,n);
D2(p2)=-1;
D2(q2)=1;

D3=zeros(1,n);
D3(p3)=-1;
D3(q3)=1;

% ţ�ٷ����N+3������
x=zeros(n+3,1);
x(1:n)=linspace(0,0.1,n);

TOL=1e-6;
k=1;

while(1)
    if k>500
        error("can not converge")
    end
    
    theta=zeros(n,1);
    F=zeros(3,1);
    
    theta=x(1:n);
    F=x(end-2:end);
    
    Rod.theta=theta;
    Rod.F=F;
    
    Rod.update;
    
    r1=Rod.K_theta*theta-transpose(Rod.Jacobian)*F;
    r2=D1*theta-Delta1;
    r3=D2*theta-Delta2;
    r4=D3*theta-Delta3;
    
    r=[r1;r2;r3;r4];
    norm(r)
    
    if(norm(r)<TOL)
        fprintf('Newton Method converge: iteration = %d\n',k-1)
        fprintf('norm(e) = %E\n',norm(r))
        break;
    end
    
    J=zeros(n+3);
    
    J(1:n,1:n)=Rod.K_theta-1*Rod.partial;
    J(1:n,n+1:end)=-transpose(Rod.Jacobian);
    J(n+1,1:n)=D1;
    J(n+2,1:n)=D2;
    J(n+3,1:n)=D3;
    
    delta=-pinv(J)*r;
    
    x=x+delta;

    k=k+1;
end

Rod.plot_all