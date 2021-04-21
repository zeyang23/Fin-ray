%  使用cosserat模型求解无外力、无刚性约束的Finray

% 21-04-21
% 求解时所用的算法很重要
% 本例中使用levenberg-marquardt无法求解，但是使用trust-region或trust-region-dogleg就可以求解


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


f=@(x) check_balance(x,psi,pA,pB,LA,LB,EA,EB,IA,IB);


x0=zeros(9,1);

options_A = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
options_B = optimoptions('fsolve','Display','off','Algorithm','trust-region');
options_C = optimoptions('fsolve','Display','off','Algorithm','trust-region-dogleg');
[x_cosserat,fval_cosserat,exitflag_cosserat,output_cosserat]=fsolve(f,x0,options_C);


%% 验证结果

na_0=x_cosserat(1:2);
ma_0=x_cosserat(3);

nb_0=x_cosserat(4:5);
mb_0=x_cosserat(6);

span_a = [0 LA];
ya_0 = [pA;na_0;ma_0];
options_a=odeset('MaxStep',1e-2);

[~,YA] = ode45(@(s,y) get_ydot(s,y,LA,EA,IA), span_a, ya_0,options_a);
ye_a=transpose(YA(end,:));


span_b = [0 LB];
yb_0 = [pB;nb_0;mb_0];
options_b=odeset('MaxStep',1e-2);
[~,YB] = ode45(@(s,y) get_ydot(s,y,LB,EB,IB), span_b, yb_0,options_b);
ye_b=transpose(YB(end,:));


hold on
plot(YA(:,1),YA(:,2))
plot(YB(:,1),YB(:,2))
axis equal


%%
function res=check_balance(x,psi,pA,pB,LA,LB,EA,EB,IA,IB)
    
    na_0=x(1:2);
    ma_0=x(3);
    
    nb_0=x(4:5);
    mb_0=x(6);
    
    pe=x(7:9);
    
    
    span_a = [0 LA];
    ya_0 = [pA;na_0;ma_0];
    options_a=odeset('MaxStep',1e-2);
    
    [~,YA] = ode45(@(s,y) get_ydot(s,y,LA,EA,IA), span_a, ya_0,options_a);
    ye_a=transpose(YA(end,:));
    
    
    span_b = [0 LB];
    yb_0 = [pB;nb_0;mb_0];
    options_b=odeset('MaxStep',1e-2);
    [~,YB] = ode45(@(s,y) get_ydot(s,y,LB,EB,IB), span_b, yb_0,options_b);
    ye_b=transpose(YB(end,:));
    
    
    res=zeros(6,1);
    res(1:2)=ye_a(1:2)-pe(1:2);
    res(3)=ye_a(3)-pe(3);
    
    res(4:5)=ye_b(1:2)-pe(1:2);
    res(6)=ye_b(3)-pe(3)-psi;
    
    res(7:8)=ye_a(4:5)+ye_b(4:5);
    res(9)=transpose(ye_a(1:2))*[0 1;-1 0]*ye_a(4:5)+ye_a(6)...
          +transpose(ye_b(1:2))*[0 1;-1 0]*ye_b(4:5)+ye_b(6);

end


function ydot=get_ydot(s,y,L,E,I)
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