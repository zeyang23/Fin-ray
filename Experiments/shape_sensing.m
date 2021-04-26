function [Rod_pos,sensorA_pos,sensorB_pos,sensorC_pos]=shape_sensing(U1,U2,U3)
wid=14e-3;
thi=0.3e-3;
E=197*1e9;

L0=0.25;
n=60;

% �������ĳ���
sensor_length=38e-3;
sensor_length_ratio=sensor_length/L0;

% ������������λ��

% �ӹ̶��˵���ˣ���������ΪA B C

sensorA_center_ratio=59e-3/L0;
sensorB_center_ratio=125e-3/L0;
sensorC_center_ratio=191e-3/L0;


pA=fix((sensorA_center_ratio-sensor_length_ratio/2)*n);
qA=fix((sensorA_center_ratio+sensor_length_ratio/2)*n);

pB=fix((sensorB_center_ratio-sensor_length_ratio/2)*n);
qB=fix((sensorB_center_ratio+sensor_length_ratio/2)*n);

pC=fix((sensorC_center_ratio-sensor_length_ratio/2)*n);
qC=fix((sensorC_center_ratio+sensor_length_ratio/2)*n);



DA=zeros(1,n);
for i=pA:qA
    DA(i)=1;
end

DB=zeros(1,n);
for i=pB:qB
    DB(i)=1;
end

DC=zeros(1,n);
for i=pC:qC
    DC(i)=1;
end


%% ʹ��fsolve���

x0=zeros(n+3,1);
options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');


coff_sensor_1=[0.876476500996772,0.050544742835559];
coff_sensor_2=[0.780420904457750,0.064426949451919];
coff_sensor_3=[1.137622996851578,0.037800943806758];


% coff_sensor_1(2)=0;
% coff_sensor_2(2)=0;
% coff_sensor_3(2)=0;


% ��λʱ�ĵ�ѹ
U1_0=1.93;
U2_0=2.79;
U3_0=1.40;

U1_delta=U1-U1_0;
U2_delta=U2-U2_0;
U3_delta=U3-U3_0;

DeltaA=polyval(coff_sensor_3,U3_delta);
DeltaB=polyval(coff_sensor_1,U1_delta);
DeltaC=polyval(coff_sensor_2,U2_delta);


f_sensor=@(x) myfun(x,E,L0,wid,thi,n,DA,DB,DC,DeltaA,DeltaB,DeltaC);
x_sensor= fsolve(f_sensor,x0,options);


Rod_solve=planar_nR(E,L0,wid,thi,n,[0;0;0]);
Rod_solve.theta=x_sensor(1:n);
Rod_solve.F=x_sensor(end-2:end);

Rod_solve.update;

Rod_solve.cal_posall;

Rod_pos=Rod_solve.pos_all;
sensorA_pos=Rod_solve.pos_all(1+pA:1+qA,:);
sensorB_pos=Rod_solve.pos_all(1+pB:1+qB,:);
sensorC_pos=Rod_solve.pos_all(1+pC:1+qC,:);

% hold on
% plot(sensor1_pos(:,1),sensor1_pos(:,2),'o')
% plot(sensor2_pos(:,1),sensor2_pos(:,2),'o')
% plot(sensor3_pos(:,1),sensor3_pos(:,2),'o')
% 
% axis([0 L0 0 3/5*L0]);
% axis equal
% 
% hold on

% Delta1_degree=DeltaA/pi*180;
% Delta2_degree=DeltaB/pi*180;
% Delta3_degree=DeltaC/pi*180;
% 
% disp("��1���н�")
% disp(Delta1_degree)
% disp("��2���н�")
% disp(Delta2_degree)
% disp("��3���н�")
% disp(Delta3_degree)
end

function [f,grad] = myfun(x,E,L0,wid,thi,n,D1,D2,D3,Delta1,Delta2,Delta3)
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
    
    f=[r1;r2;r3;r4];
    
    if nargout > 1
        J=zeros(n+3);
    
        J(1:n,1:n)=Rod.K_theta-1*Rod.partial;
        J(1:n,n+1:end)=-transpose(Rod.Jacobian);
        J(n+1,1:n)=D1;
        J(n+2,1:n)=D2;
        J(n+3,1:n)=D3;
        grad=J;
    end
end