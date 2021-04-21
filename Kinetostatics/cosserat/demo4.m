% 使用cosserat模型求解无外力、有1根刚性约束的Finray

% 求解失败，狄拉克函数处理的估计不太对


%% Cosserat Model shooting method
clear
clc

wid=14e-3;
thi=0.3e-3;
E=197*1e9;

Iz=1/12*thi.^3*wid;

pdes=[0;0;pi];

L0=1;

LA=0.8*L0;
LB=L0;

EA=E;
EB=E;

IA=Iz;
IB=Iz;


xA=0;
yA=0;

xB=0.35*L0;
yB=0*L0;


psi_degree=20;
alpha_degree=80;
beta_degree=100;

psi=psi_degree/180*pi;
alpha=alpha_degree/180*pi;
beta=beta_degree/180*pi;

pA=[xA;yA;alpha];
pB=[xB;yB;beta];

Lcon=1/2*sqrt((xA-xB)^2+(yA-yB)^2);
lambdaA=1/2*LA;
lambdaB=1/2*LB;

f=@(x) check_balance(x,Lcon,lambdaA,lambdaB,psi,pA,pB,LA,LB,EA,EB,IA,IB);


x0=zeros(11,1);

options_A = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
options_B = optimoptions('fsolve','Display','off','Algorithm','trust-region');
options_C = optimoptions('fsolve','Display','off','Algorithm','trust-region-dogleg');
[X_cosserat,fval_cosserat,exitflag_cosserat,output_cosserat]=fsolve(f,x0,options_B);


%% 验证结果

na_0=X_cosserat(1:2);
ma_0=X_cosserat(3);

nb_0=X_cosserat(4:5);
mb_0=X_cosserat(6);

fcon=X_cosserat(10);
gamma=X_cosserat(11);

span_a = [0 LA];
ya_0 = [pA;na_0;ma_0];
options_a=odeset('MaxStep',1e-2);
FA=[-fcon*cos(gamma);-fcon*sin(gamma)];
[~,YA] = ode45(@(s,y) get_ydot(s,y,lambdaA,FA,LA,EA,IA), span_a, ya_0,options_a);
ye_a=transpose(YA(end,:));


span_b = [0 LB];
yb_0 = [pB;nb_0;mb_0];
options_b=odeset('MaxStep',1e-2);
FB=-FA;
[~,YB] = ode45(@(s,y) get_ydot(s,y,lambdaB,FB,LB,EB,IB), span_b, yb_0,options_b);
ye_b=transpose(YB(end,:));


hold on
plot(YA(:,1),YA(:,2))
plot(YB(:,1),YB(:,2))
axis equal


%%
function res=check_balance(x,Lcon,lambdaA,lambdaB,psi,pA,pB,LA,LB,EA,EB,IA,IB)
    
    na_0=x(1:2);
    ma_0=x(3);
    
    nb_0=x(4:5);
    mb_0=x(6);
    
    pe=x(7:9);
    
    fcon=x(10);
    gamma=x(11);
    
    
    span_a = [0 LA];
    ya_0 = [pA;na_0;ma_0];
    options_a=odeset('MaxStep',1e-2);
    FA=[-fcon*cos(gamma);-fcon*sin(gamma)];
    [~,YA] = ode45(@(s,y) get_ydot(s,y,lambdaA,FA,LA,EA,IA), span_a, ya_0,options_a);
    ye_a=transpose(YA(end,:));
    
    
    span_b = [0 LB];
    yb_0 = [pB;nb_0;mb_0];
    options_b=odeset('MaxStep',1e-2);
    FB=-FA;
    [~,YB] = ode45(@(s,y) get_ydot(s,y,lambdaB,FB,LB,EB,IB), span_b, yb_0,options_b);
    ye_b=transpose(YB(end,:));
    
    
    res=zeros(11,1);
    res(1:2)=ye_a(1:2)-pe(1:2);
    res(3)=ye_a(3)-pe(3);
    
    res(4:5)=ye_b(1:2)-pe(1:2);
    res(6)=ye_b(3)-pe(3)-psi;
    
    res(7:8)=ye_a(4:5)+ye_b(4:5);
    res(9)=transpose(ye_a(1:2))*[0 1;-1 0]*ye_a(4:5)+ye_a(6)...
          +transpose(ye_b(1:2))*[0 1;-1 0]*ye_b(4:5)+ye_b(6);
      
    res(10:11)=ye_a(1:2)+[Lcon*cos(gamma);Lcon*sin(gamma)]-ye_b(1:2);

end


function ydot=get_ydot(s,y,lambda,F,L,E,I)
    ydot=zeros(size(y));
    
    r=y(1:2);
    theta=y(3);
    n=y(4:5);
    m=y(6);
    
    delta=1e-2;
    if norm(s-lambda)<=delta/2
        f=F/delta; 
    else
        f=[0;0];
    end
    
    l=0;

    ydot(1:2)=[cos(theta);sin(theta)];
    ydot(3)=1/(E*I)*m;
    ydot(4:5)=-f;
    ydot(6)=-[cos(theta) sin(theta)]*[0 1;-1 0]*n-l;
end