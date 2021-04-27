% 基于cosserat model
% 单个柔性板的shape sensing
% 考虑传感器的模型
% 认为液态金属传感器的阻值仅与拉伸有关，不受弯曲影响


%% 给定末端位姿，正向计算
clear
clc

L0=0.25;

% pdes=[L0;0;0];
pdes=[0.6*L0;0.6*L0;pi/3];


% 3个传感器的位置
SA_start=0.040;
SA_end=0.078;

SB_start=0.106;
SB_end=0.144;

SC_start=0.172;
SC_end=0.210;

wid=14e-3;
thi=0.3e-3;
E=197*1e9;

Iz=1/12*thi.^3*wid;


f=@(x) check_balance(x,L0,E,Iz,pdes);


x0=zeros(6,1);

options_a = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
options_b = optimoptions('fsolve','Display','off','Algorithm','trust-region');
options_c = optimoptions('fsolve','Display','off','Algorithm','trust-region-dogleg');
[x_cosserat,fval_cosserat,exitflag_cosserat,output_cosserat]=fsolve(f,x0,options_a);


% 验证结果

span = [0 L0];
y0 = [0;0;0;x_cosserat(1:3)];
options3=odeset('MaxStep',1e-2);
sol_Y = ode45(@(s,y) get_ydot(s,y,L0,E,Iz,x_cosserat(4:6)), span, y0,options3);

hold on
plot(sol_Y.y(1,:),sol_Y.y(2,:),'r')


N=5;
H=0.012;


lengthA=0;
SA_series=linspace(SA_start,SA_end,N);
posA_rod=zeros(N,2);
posA_sensor=zeros(N,2);
for i=1:N
    YA_now=deval(sol_Y,SA_series(i));
    posA_rod(i,:)=[YA_now(1),YA_now(2)];
    posA_sensor(i,:)=posA_rod(i,:)+[H*sin(YA_now(3)),-H*cos(YA_now(3))];
end
for i=2:N
    lengthA=lengthA+norm(posA_sensor(i,:)-posA_sensor(i-1,:));
end

lengthB=0;
SB_series=linspace(SB_start,SB_end,N);
posB_rod=zeros(N,2);
posB_sensor=zeros(N,2);
for i=1:N
    YB_now=deval(sol_Y,SB_series(i));
    posB_rod(i,:)=[YB_now(1),YB_now(2)];
    posB_sensor(i,:)=posB_rod(i,:)+[H*sin(YB_now(3)),-H*cos(YB_now(3))];
end
for i=2:N
    lengthB=lengthB+norm(posB_sensor(i,:)-posB_sensor(i-1,:));
end

lengthC=0;
SC_series=linspace(SC_start,SC_end,N);
posC_rod=zeros(N,2);
posC_sensor=zeros(N,2);
for i=1:N
    YC_now=deval(sol_Y,SC_series(i));
    posC_rod(i,:)=[YC_now(1),YC_now(2)];
    posC_sensor(i,:)=posC_rod(i,:)+[H*sin(YC_now(3)),-H*cos(YC_now(3))];
end
for i=2:N
    lengthC=lengthC+norm(posC_sensor(i,:)-posC_sensor(i-1,:));
end





%% 根据曲率反解柔性板形状
g=@(x) check_shape(x,L0,E,Iz,lengthA,lengthB,lengthC,SA_start,SA_end,SB_start,SB_end,SC_start,SC_end);


x0=zeros(6,1);

options_a = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
options_b = optimoptions('fsolve','Display','off','Algorithm','trust-region');
options_c = optimoptions('fsolve','Display','off','Algorithm','trust-region-dogleg');
[x_shape,fval_shape,exitflag_shape,output_shape]=fsolve(g,x0,options_a);

span = [0 L0];
y0 = [0;0;0;x_shape(1:3)];
options3=odeset('MaxStep',1e-2);
sol_Y2 = ode45(@(s,y) get_ydot(s,y,L0,E,Iz,x_shape(4:6)), span, y0,options3);


figure

hold on
plot(sol_Y2.y(1,:),sol_Y2.y(2,:))


N=5;
H=0.012;


lengthA=0;
SA_series=linspace(SA_start,SA_end,N);
posA_rod=zeros(N,2);
posA_sensor=zeros(N,2);
for i=1:N
    YA_now=deval(sol_Y2,SA_series(i));
    posA_rod(i,:)=[YA_now(1),YA_now(2)];
    posA_sensor(i,:)=posA_rod(i,:)+[H*sin(YA_now(3)),-H*cos(YA_now(3))];
    
    plot([posA_rod(i,1) posA_sensor(i,1)],[posA_rod(i,2) posA_sensor(i,2)])
end
plot(posA_sensor(:,1),posA_sensor(:,2),'-o','MarkerSize',2)
for i=2:N
    lengthA=lengthA+norm(posA_sensor(i,:)-posA_sensor(i-1,:));
end

lengthB=0;
SB_series=linspace(SB_start,SB_end,N);
posB_rod=zeros(N,2);
posB_sensor=zeros(N,2);
for i=1:N
    YB_now=deval(sol_Y2,SB_series(i));
    posB_rod(i,:)=[YB_now(1),YB_now(2)];
    posB_sensor(i,:)=posB_rod(i,:)+[H*sin(YB_now(3)),-H*cos(YB_now(3))];
    plot([posB_rod(i,1) posB_sensor(i,1)],[posB_rod(i,2) posB_sensor(i,2)])
end
plot(posB_sensor(:,1),posB_sensor(:,2),'-o','MarkerSize',2)
for i=2:N
    lengthB=lengthB+norm(posB_sensor(i,:)-posB_sensor(i-1,:));
end

lengthC=0;
SC_series=linspace(SC_start,SC_end,N);
posC_rod=zeros(N,2);
posC_sensor=zeros(N,2);
for i=1:N
    YC_now=deval(sol_Y2,SC_series(i));
    posC_rod(i,:)=[YC_now(1),YC_now(2)];
    posC_sensor(i,:)=posC_rod(i,:)+[H*sin(YC_now(3)),-H*cos(YC_now(3))];
    plot([posC_rod(i,1) posC_sensor(i,1)],[posC_rod(i,2) posC_sensor(i,2)])
end
plot(posC_sensor(:,1),posC_sensor(:,2),'-o','MarkerSize',2)
for i=2:N
    lengthC=lengthC+norm(posC_sensor(i,:)-posC_sensor(i-1,:));
end

axis equal
grid on

%%
function res=check_shape(x,L,E,I,lengthA,lengthB,lengthC,SA_start,SA_end,SB_start,SB_end,SC_start,SC_end)
    n0=x(1:2);
    m0=x(3);
    Fe=x(4:6);
    
    span = [0 L];
    y0 = [0;0;0;n0;m0];
    
    options=odeset('MaxStep',1e-2);
    
    sol_Y = ode45(@(s,y) get_ydot(s,y,L,E,I,Fe), span, y0,options);
    
    ye=sol_Y.y(:,end);
    
    
    N=5;
    H=0.012;


    lengthA_real=0;
    SA_series=linspace(SA_start,SA_end,N);
    posA_rod=zeros(N,2);
    posA_sensor=zeros(N,2);
    for i=1:N
        YA_now=deval(sol_Y,SA_series(i));
        posA_rod(i,:)=[YA_now(1),YA_now(2)];
        posA_sensor(i,:)=posA_rod(i,:)+[H*sin(YA_now(3)),-H*cos(YA_now(3))];
    end
    for i=2:N
        lengthA_real=lengthA_real+norm(posA_sensor(i,:)-posA_sensor(i-1,:));
    end

    lengthB_real=0;
    SB_series=linspace(SB_start,SB_end,N);
    posB_rod=zeros(N,2);
    posB_sensor=zeros(N,2);
    for i=1:N
        YB_now=deval(sol_Y,SB_series(i));
        posB_rod(i,:)=[YB_now(1),YB_now(2)];
        posB_sensor(i,:)=posB_rod(i,:)+[H*sin(YB_now(3)),-H*cos(YB_now(3))];
    end
    for i=2:N
        lengthB_real=lengthB_real+norm(posB_sensor(i,:)-posB_sensor(i-1,:));
    end

    lengthC_real=0;
    SC_series=linspace(SC_start,SC_end,N);
    posC_rod=zeros(N,2);
    posC_sensor=zeros(N,2);
    for i=1:N
        YC_now=deval(sol_Y,SC_series(i));
        posC_rod(i,:)=[YC_now(1),YC_now(2)];
        posC_sensor(i,:)=posC_rod(i,:)+[H*sin(YC_now(3)),-H*cos(YC_now(3))];
    end
    for i=2:N
        lengthC_real=lengthC_real+norm(posC_sensor(i,:)-posC_sensor(i-1,:));
    end

    
    res=zeros(6,1);
    res(1)=lengthA-lengthA_real;
    res(2)=lengthB-lengthB_real;
    res(3)=lengthC-lengthC_real;
    res(4:5)=ye(4:5)-Fe(1:2);
    res(6)=transpose(ye(1:2))*[0 1;-1 0]*ye(4:5)+ye(6)-transpose(ye(1:2))*[0 1;-1 0]*Fe(1:2)-Fe(3);

end

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