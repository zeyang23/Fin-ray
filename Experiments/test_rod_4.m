% 生成标定所需的末端位姿参数

clear
clc

wid=14e-3;
thi=0.3e-3;
E=197*1e9;

L0=0.25;
n=100;

pdes1=[0.45*L0;0.45*L0;pi];
pdes2=[0.45*L0;0.45*L0;5*pi/6];

pdes3=[0.9*L0;0.05*L0;pi/6];
pdes4=[0.9*L0;0.05*L0;pi/3];

pdes5=[0*L0;0.4*L0;pi];
pdes6=[0*L0;0.8*L0;pi];

pdes7=[0.1*L0;0.4*L0;pi];
pdes8=[0.1*L0;0.8*L0;pi];

pdes9=[0.2*L0;0.4*L0;pi];
pdes10=[0.2*L0;0.8*L0;pi];

R1=planar_nR(E,L0,wid,thi,n,pdes2);


f_rod=@(x) cal_balance(x,R1);
options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');
x0=zeros(R1.n_seg+3,1);
[x_solve_rod,fval_rod,exitflag_rod,output_rod] = fsolve(f_rod,x0,options);

R1.theta=x_solve_rod(1:R1.n_seg);
R1.F=x_solve_rod(R1.n_seg+1:R1.n_seg+3);

R1.update;
R1.plot_all;


% 传感器

% 传感器的长度
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

Delta1=sum(R1.theta(p1:q1));
Delta2=sum(R1.theta(p2:q2));
Delta3=sum(R1.theta(p3:q3));

Delta1_degree=Delta1/pi*180;
Delta2_degree=Delta2/pi*180;
Delta3_degree=Delta3/pi*180;

disp("第1个夹角")
disp(Delta1_degree)
disp("第2个夹角")
disp(Delta2_degree)
disp("第3个夹角")
disp(Delta3_degree)

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

% 画出3个传感器的位置
R1.cal_posall;
sensor1_pos=R1.pos_all(1+p1:1+q1,:);
sensor2_pos=R1.pos_all(1+p2:1+q2,:);
sensor3_pos=R1.pos_all(1+p3:1+q3,:);

hold on
plot(sensor1_pos(:,1),sensor1_pos(:,2),'o')
plot(sensor2_pos(:,1),sensor2_pos(:,2),'o')
plot(sensor3_pos(:,1),sensor3_pos(:,2),'o')


function [r,J]=cal_balance(x,Rod)
    Rod.theta=x(1:Rod.n_seg);
    Rod.F=x(Rod.n_seg+1:Rod.n_seg+3);
    
    Rod.update;
    r=Rod.cal_r;
    J=Rod.cal_rdot;
end