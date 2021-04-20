% 对比主轴分解法与Cosserat法

% 21-04-18
% 在一些情况下，两种方法求得的形状是相同的，即使形状相同，求得的末端的力的大小也不同，令人迷惑
% 在另一些情况下，两种方法求得的形状不同

%% 主轴分解

clear
clc

wid=14e-3;
thi=0.3e-3;
E=197*1e9;

L0=1;
n=50;

pdes=[0.6*L0;0.6*L0;pi/2];

R1=planar_nR(E,L0,wid,thi,n,pdes);


f=@(x) cal_balance(x,R1);
options1 = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');
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

options2 = optimoptions('fsolve','Display','off');
x_cosserat=fsolve(f,x0,options2);


% 验证结果

span = [0 L];
y0 = [0;0;0;x_cosserat(1:3)];
options3=odeset('MaxStep',1e-2);
[s,Y] = ode23s(@(s,y) get_ydot(s,y,L,E,Iz,x_cosserat(4:6)), span, y0,options3);

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
    
    [~,Y] = ode23s(@(s,y) get_ydot(s,y,L,E,I,Fe), span, y0,options);
    
    ye=transpose(Y(end,:));
    
    res=zeros(6,1);
    res(1:2)=ye(1:2)-pdes(1:2);
    res(3)=ye(3)-pdes(3);
    res(4:5)=ye(4:5)-Fe(1:2);
    res(6)=-[cos(ye(3)) sin(ye(3))]*[0 1;-1 0]*ye(4:5)+ye(6)-transpose(ye(1:2))*[0 1;-1 0]*Fe(1:2)-Fe(3);
    
end


function ydot=get_ydot(s,y,L,E,I,Fe)
    ydot=zeros(size(y));
    
    r=y(1:2);
    theta=y(3);
    n=y(4:5);
    m=y(6);
    
    delta=1e-16;
    if (L-s)<delta
        f=Fe(1:2);
        l=Fe(3);
    else
        f=[0;0];
        l=0;
    end
    
    ydot(1:2)=[cos(theta);sin(theta)];
    ydot(3)=1/(E*I)*m;
    ydot(4:5)=f;
    ydot(6)=-[cos(theta) sin(theta)]*[0 1;-1 0]*n-l;
end

function [r,J]=cal_balance(x,Rod)
    Rod.theta=x(1:Rod.n_seg);
    Rod.F=x(Rod.n_seg+1:Rod.n_seg+3);
    
    Rod.update;
    r=Rod.cal_r;
    J=Rod.cal_rdot;
end