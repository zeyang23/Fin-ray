% 近似为圆，给圆的半径，再转换成角度

clear
clc

wid=14e-3;
thi=0.3e-3;
E=197*1e9;

L0=0.25;
n=100;

sensor_length=38e-3;
sensor_length_ratio=sensor_length/L0;

% 传感器的中心位置
sensor1_center_ratio=1/4;
sensor2_center_ratio=2/4;
sensor3_center_ratio=3/4;


p1=fix((sensor1_center_ratio-sensor_length_ratio/2)*n);
q1=fix((sensor1_center_ratio+sensor_length_ratio/2)*n);

p2=fix((sensor2_center_ratio-sensor_length_ratio/2)*n);
q2=fix((sensor2_center_ratio+sensor_length_ratio/2)*n);

p3=fix((sensor3_center_ratio-sensor_length_ratio/2)*n);
q3=fix((sensor3_center_ratio+sensor_length_ratio/2)*n);

% Delta1=sum(R1.theta(p1:q1));
% Delta2=sum(R1.theta(p2:q2));
% Delta3=sum(R1.theta(p3:q3));
% 
% Delta1_degree=Delta1/pi*180;
% Delta2_degree=Delta2/pi*180;
% Delta3_degree=Delta3/pi*180;
% 
% disp("第1个夹角")
% disp(Delta1_degree)
% disp("第2个夹角")
% disp(Delta2_degree)
% disp("第3个夹角")
% disp(Delta3_degree)

D1=zeros(1,n);
for i=p1:q1
    D1(i)=1;
end

D2=zeros(1,n);
for i=p2:q2
    D2(i)=1;
end

D3=zeros(1,n);
for i=p3:q3
    D3(i)=1;
end


%% 使用fsolve求解

x0=zeros(n+3,1);
options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

L=0.038;

Delta1=L/(0.5/pi);
Delta2=L/(0.5/pi);
Delta3=L/(0.5/pi);


f_sensor=@(x) myfun(x,E,L0,wid,thi,n,D1,D2,D3,Delta1,Delta2,Delta3);
[x_sensor,fval_sensor,exitflag_sensor,output_sensor] = fsolve(f_sensor,x0,options);


Rod_solve=planar_nR(E,L0,wid,thi,n,[0;0;0]);
Rod_solve.theta=x_sensor(1:n);
Rod_solve.F=x_sensor(end-2:end);

Rod_solve.update;
Rod_solve.plot_all;

Rod_solve.cal_posall;
sensor1_pos=Rod_solve.pos_all(1+p1:1+q1,:);
sensor2_pos=Rod_solve.pos_all(1+p2:1+q2,:);
sensor3_pos=Rod_solve.pos_all(1+p3:1+q3,:);

hold on
plot(sensor1_pos(:,1),sensor1_pos(:,2),'o')
plot(sensor2_pos(:,1),sensor2_pos(:,2),'o')
plot(sensor3_pos(:,1),sensor3_pos(:,2),'o')

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