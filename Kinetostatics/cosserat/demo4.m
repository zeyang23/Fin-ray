% 使用cosserat模型求解无外力、有1根刚性约束的Finray
% 尝试使用分段的方法来处理

% 21-04-21
% 分段处理的效果不错，能够求解

%% Cosserat Model shooting method
clear
clc

wid=14e-3;
thi=0.3e-3;
E=197*1e9;

Iz=1/12*thi.^3*wid;


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
[x_cosserat,fval_cosserat,exitflag_cosserat,output_cosserat]=fsolve(f,x0,options_B);



%% 验证求解结果
na_0=x_cosserat(1:2);
ma_0=x_cosserat(3);

nb_0=x_cosserat(4:5);
mb_0=x_cosserat(6);

pe=x_cosserat(7:9);

fcon=x_cosserat(10);
gamma=x_cosserat(11);
    
    
span_a1 = [0 lambdaA];
ya1_0 = [pA;na_0;ma_0];
options_a=odeset('MaxStep',1e-2);
[~,YA_in] = ode45(@(s,y) get_ydot(s,y,LA,EA,IA), span_a1, ya1_0,options_a);
yin_a=transpose(YA_in(end,:));

FA=[-fcon*cos(gamma);-fcon*sin(gamma)];
span_a2 = [lambdaA,LA];
ya2_0=yin_a+[0;0;0;-FA;0];
[~,YA] = ode45(@(s,y) get_ydot(s,y,LA,EA,IA), span_a2, ya2_0,options_a);
ye_a=transpose(YA(end,:));


span_b1 = [0 lambdaB];
yb1_0 = [pB;nb_0;mb_0];
options_b=odeset('MaxStep',1e-2);
[~,YB_in] = ode45(@(s,y) get_ydot(s,y,LB,EB,IB), span_b1, yb1_0,options_b);
yin_b=transpose(YB_in(end,:));

FB=-FA;
span_b2 = [lambdaB,LB];
yb2_0=yin_b+[0;0;0;-FB;0];
[~,YB] = ode45(@(s,y) get_ydot(s,y,LA,EA,IA), span_b2, yb2_0,options_b);
ye_b=transpose(YB(end,:));

YA_total=[YA_in;YA(2:end,:)];
YB_total=[YB_in;YB(2:end,:)];

hold on
plot(YA_total(:,1),YA_total(:,2))
plot(YB_total(:,1),YB_total(:,2))
axis equal

plot([yin_a(1) yin_b(1)],[yin_a(2) yin_b(2)])


%%
function res=check_balance(x,Lcon,lambdaA,lambdaB,psi,pA,pB,LA,LB,EA,EB,IA,IB)
    
    na_0=x(1:2);
    ma_0=x(3);
    
    nb_0=x(4:5);
    mb_0=x(6);
    
    pe=x(7:9);
    
    fcon=x(10);
    gamma=x(11);
    
    
    span_a1 = [0 lambdaA];
    ya1_0 = [pA;na_0;ma_0];
    options_a=odeset('MaxStep',1e-2);
    [~,YA_in] = ode45(@(s,y) get_ydot(s,y,LA,EA,IA), span_a1, ya1_0,options_a);
    yin_a=transpose(YA_in(end,:));
    
    FA=[-fcon*cos(gamma);-fcon*sin(gamma)];
    span_a2 = [lambdaA,LA];
    ya2_0=yin_a+[0;0;0;-FA;0];
    [~,YA] = ode45(@(s,y) get_ydot(s,y,LA,EA,IA), span_a2, ya2_0,options_a);
    ye_a=transpose(YA(end,:));
    
    
    span_b1 = [0 lambdaB];
    yb1_0 = [pB;nb_0;mb_0];
    options_b=odeset('MaxStep',1e-2);
    [~,YB_in] = ode45(@(s,y) get_ydot(s,y,LB,EB,IB), span_b1, yb1_0,options_b);
    yin_b=transpose(YB_in(end,:));
    
    FB=-FA;
    span_b2 = [lambdaB,LB];
    yb2_0=yin_b+[0;0;0;-FB;0];
    [~,YB] = ode45(@(s,y) get_ydot(s,y,LA,EA,IA), span_b2, yb2_0,options_b);
    ye_b=transpose(YB(end,:));
    
    
    res=zeros(11,1);
    res(1:2)=ye_a(1:2)-pe(1:2);
    res(3)=ye_a(3)-pe(3);
    
    res(4:5)=ye_b(1:2)-pe(1:2);
    res(6)=ye_b(3)-pe(3)-psi;
    
    res(7:8)=ye_a(4:5)+ye_b(4:5);
    res(9)=transpose(ye_a(1:2))*[0 1;-1 0]*ye_a(4:5)+ye_a(6)...
          +transpose(ye_b(1:2))*[0 1;-1 0]*ye_b(4:5)+ye_b(6);
      
    res(10:11)=yin_a(1:2)+[Lcon*cos(gamma);Lcon*sin(gamma)]-yin_b(1:2);

end


function ydot=get_ydot(s,y,L,E,I)
    ydot=zeros(size(y));
    
    r=y(1:2);
    theta=y(3);
    n=y(4:5);
    m=y(6);
    
    f=[0;0];
    l=0;

    ydot(1:2)=[cos(theta);sin(theta)];
    ydot(3)=1/(E*I)*m;
    ydot(4:5)=-f;
    ydot(6)=-[cos(theta) sin(theta)]*[0 1;-1 0]*n-l;
end