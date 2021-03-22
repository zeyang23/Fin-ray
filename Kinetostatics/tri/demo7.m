% ���ԣ�ͨ������3�������������ݷ������԰����״
% ʹ��fsolve�����С�������⣬ʹ�õ���ʵ��Levenberg-Marquardt algorithm����nonlsq����һ����

% 21-03-22 4����������ĳЩ�������������
% ���ṩ�ݶ��ǿ��Խ�ģ��ṩ�ݶȺ������ٶ�������졣ע�⣬�����ṩ���ݶ�Ҫ��partial�����CheckGradient�ᱨ��
% ���ڳ�ֵ�����ԡ����ڽ������ֵ���ԭʼֵ���ƫ��Ŀ����ԡ�


% ����4��Delta
clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=40;


pdes=[0.6*L0;0.6*L0;0];   % �����ĩ�˲������ǳɹ���

R1=planar_nR(E,L0,wid,thi,n,pdes);

TOL=1e-6;
R1.Newton(TOL);
R1.plot_all;
hold on

% �������Ĵ�С��λ��
sensor_length=1/8;

p1=fix((1/5-sensor_length/2)*n);
q1=fix((1/5+sensor_length/2)*n);

p2=fix((2/5-sensor_length/2)*n);
q2=fix((2/5+sensor_length/2)*n);

p3=fix((3/5-sensor_length/2)*n);
q3=fix((3/5+sensor_length/2)*n);

p4=fix((4/5-sensor_length/2)*n);
q4=fix((4/5+sensor_length/2)*n);

Delta1=R1.theta(q1)-R1.theta(p1);
Delta2=R1.theta(q2)-R1.theta(p2);
Delta3=R1.theta(q3)-R1.theta(p3);
Delta4=R1.theta(q4)-R1.theta(p4);

D1=zeros(1,n);
D1(p1)=-1;
D1(q1)=1;

D2=zeros(1,n);
D2(p2)=-1;
D2(q2)=1;

D3=zeros(1,n);
D3(p3)=-1;
D3(q3)=1;

D4=zeros(1,n);
D4(p4)=-1;
D4(q4)=1;



% ʹ��fsolve���

x0=zeros(n+3,1);
options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',true);

f=@(x) myfun(x,E,L0,wid,thi,n,D1,D2,D3,D4,Delta1,Delta2,Delta3,Delta4);
[x,fval,exitflag,output] = fsolve(f,x0,options);
% [x,fval,exitflag,output]= fsolve(f,x0);

theta_real=R1.theta;
x_real=[R1.theta;R1.F];

disp(norm(x(1:n)-theta_real));
disp(norm(x-x_real));

Rod_solve=planar_nR(E,L0,wid,thi,n,pdes);
Rod_solve.theta=x(1:n);
Rod_solve.F=x(end-2:end);

Rod_solve.update;
Rod_solve.plot_all;

function [f,grad] = myfun(x,E,L0,wid,thi,n,D1,D2,D3,D4,Delta1,Delta2,Delta3,Delta4)
    Rod=planar_nR(E,L0,wid,thi,n,[0;0;0]);
    
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
    r5=D4*theta-Delta4;
    
    f=[r1;r2;r3;r4;r5];
    
    if nargout > 1
        J=zeros(n+3);
    
        J(1:n,1:n)=Rod.K_theta-1*Rod.partial;
        J(1:n,n+1:end)=-transpose(Rod.Jacobian);
        J(n+1,1:n)=D1;
        J(n+2,1:n)=D2;
        J(n+3,1:n)=D3;
        J(n+4,1:n)=D4;
        grad=J;
    end
end