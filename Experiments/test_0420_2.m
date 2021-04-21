clear
clc

%% 录入实验数据

A=zeros(15,6);
A(1,:)=[59.9 65.4 157.9 2.28 2.30 2.40];
A(2,:)=[28.5 58.2 182.3 2.34 2.35 2.41];
A(3,:)=[-9.1 44.5 196.5 2.36 2.39 2.44];
A(4,:)=[-53.5 20 223.8 2.37 2.40 2.47];
A(5,:)=[-53.5 20 240.5 2.40 2.34 2.44];
A(6,:)=[-53.5 6.5 236.2 2.44 2.46 2.47];
A(7,:)=[-53.5 6.5 221.9 2.38 2.51 2.50];
A(8,:)=[-57 -12.6 256.5 2.51 2.48 2.47];
A(9,:)=[-57 -12.6 231.1 2.41 2.58 2.51];
A(10,:)=[-73.5 -22.5 231.1 2.36 2.57 2.55];
A(11,:)=[-73.5 -22.5 255.5 2.46 2.48 2.51];
A(12,:)=[-82.5 -42.5 273 2.52 2.51 2.50];
A(13,:)=[-82.5 -62 280 2.55 2.55 2.51];
A(14,:)=[-82.5 -62 308 2.61 2.43 2.48];
A(15,:)=[-82.5 -62 342 2.53 2.35 2.46];

temp=A(:,1);
A(:,1)=A(:,2);
A(:,2)=temp;

% 减掉0位的值
B=zeros(size(A));
for i=1:size(A,1)
    B(i,:)=A(i,:)-A(1,:);
end

B(:,1)=B(:,1);
B(:,2)=-B(:,2);
B(:,3)=B(:,3)/180*pi;


pdes_series=B(:,1:3);
R_series=B(:,4:end);

R_series=fliplr(R_series);


%% 计算每个位姿下的三个角度

wid=14e-3;
thi=0.3e-3;
E=197*1e9;

L0=0.25;
n=50;

pdes_series(:,1)=pdes_series(:,1)/1000+L0;
pdes_series(:,2)=pdes_series(:,2)/1000;

Delta_series=zeros(size(R_series));

for i=1:size(pdes_series,1)

    pdes=transpose(pdes_series(i,:));

    R1=planar_nR(E,L0,wid,thi,n,pdes);


    f_rod=@(x) cal_balance(x,R1);
    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off','Algorithm','levenberg-marquardt');
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
    sensor1_center_ratio=59e-3/L0;
    sensor2_center_ratio=125e-3/L0;
    sensor3_center_ratio=191e-3/L0;


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
    
    Delta_series(i,1)=Delta1_degree;
    Delta_series(i,2)=Delta2_degree;
    Delta_series(i,3)=Delta3_degree;
    
    % 画出3个传感器的位置
    R1.cal_posall;
    sensor1_pos=R1.pos_all(1+p1:1+q1,:);
    sensor2_pos=R1.pos_all(1+p2:1+q2,:);
    sensor3_pos=R1.pos_all(1+p3:1+q3,:);

    hold on
    plot(sensor1_pos(:,1),sensor1_pos(:,2),'o')
    plot(sensor2_pos(:,1),sensor2_pos(:,2),'o')
    plot(sensor3_pos(:,1),sensor3_pos(:,2),'o')
    
    axis([0 L0 -0.3*L0 0.7*L0])
    axis equal
    
    pause(0.02)
    
    if i~=size(pdes_series,1)
        clf;
    end
end


%% 画出角度与电阻的关系

figure
subplot(3,1,1)
R1_series=R_series(:,1);
Delta1_series=Delta_series(:,1);
RandD_1=sortrows([Delta1_series,R1_series],1);
plot(RandD_1(:,1),RandD_1(:,2),'-o','MarkerSize',2)
xlabel('delta')
ylabel('R')


subplot(3,1,2)
R2_series=R_series(:,2);
Delta2_series=Delta_series(:,2);
RandD_2=sortrows([Delta2_series,R2_series],1);
plot(RandD_2(:,1),RandD_2(:,2),'-o','MarkerSize',2)
xlabel('delta')
ylabel('R')


subplot(3,1,3)
R3_series=R_series(:,3);
Delta3_series=Delta_series(:,3);
RandD_3=sortrows([Delta3_series,R3_series],1);
plot(RandD_3(:,1),RandD_3(:,2),'-o','MarkerSize',2)
xlabel('delta')
ylabel('R')


%%
function [r,J]=cal_balance(x,Rod)
    Rod.theta=x(1:Rod.n_seg);
    Rod.F=x(Rod.n_seg+1:Rod.n_seg+3);
    
    Rod.update;
    r=Rod.cal_r;
    J=Rod.cal_rdot;
end