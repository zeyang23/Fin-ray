% ʹ��cosseratģ���������������1������Լ����Finray
% ����Ӵ���ΪһС�ε��߽Ӵ�

% 21-04-22
% ʹ�þ���������dirac����������⣬����³���ԱȽϲ�
% ���Ҿ��εĿ�ȸ���ôȡҲ�Ǹ�����


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
[x_cosserat,fval_cosserat,exitflag_cosserat,output_cosserat]=fsolve(f,x0,options_C);



%% ��֤�����
na_0=x_cosserat(1:2);
ma_0=x_cosserat(3);

nb_0=x_cosserat(4:5);
mb_0=x_cosserat(6);

pe=x_cosserat(7:9);

fcon=x_cosserat(10);
gamma=x_cosserat(11);
    
    
span_a = [0 LA];
ya_0 = [pA;na_0;ma_0];
options_a=odeset('MaxStep',1e-2);

FA=[-fcon*cos(gamma);-fcon*sin(gamma)];

sol_YA = ode45(@(s,y) get_ydot(s,y,lambdaA,FA,LA,EA,IA), span_a, ya_0,options_a);


span_b = [0 LB];
yb_0 = [pB;nb_0;mb_0];
options_b=odeset('MaxStep',1e-2);

FB=-FA;

sol_YB = ode45(@(s,y) get_ydot(s,y,lambdaB,FB,LB,EB,IB), span_b, yb_0,options_b);
    
yin_a=deval(sol_YA,lambdaA);
yin_b=deval(sol_YB,lambdaB);

hold on
plot(sol_YA.y(1,:),sol_YA.y(2,:))
plot(sol_YB.y(1,:),sol_YB.y(2,:))
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
    
    
    span_a = [0 LA];
    ya_0 = [pA;na_0;ma_0];
    options_a=odeset('MaxStep',1e-2);

    FA=[-fcon*cos(gamma);-fcon*sin(gamma)];
    
    sol_YA = ode45(@(s,y) get_ydot(s,y,lambdaA,FA,LA,EA,IA), span_a, ya_0,options_a);
    
    
    span_b = [0 LB];
    yb_0 = [pB;nb_0;mb_0];
    options_b=odeset('MaxStep',1e-2);
    
    FB=-FA;
    
    sol_YB = ode45(@(s,y) get_ydot(s,y,lambdaB,FB,LB,EB,IB), span_b, yb_0,options_b);
    
    
    ye_a=sol_YA.y(:,end);
    ye_b=sol_YB.y(:,end);
    
    yin_a=deval(sol_YA,lambdaA);
    yin_b=deval(sol_YB,lambdaB);
    
    
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


function ydot=get_ydot(s,y,lambda,F,L,E,I)
    ydot=zeros(size(y));
    
    r=y(1:2);
    theta=y(3);
    n=y(4:5);
    m=y(6);
    
    l=0;
    
    delta=1e-2;
    if norm(s-lambda)<=delta/2
        f=F/delta;
    else
        f=[0;0];
    end

    ydot(1:2)=[cos(theta);sin(theta)];
    ydot(3)=1/(E*I)*m;
    ydot(4:5)=-f;
    ydot(6)=-[cos(theta) sin(theta)]*[0 1;-1 0]*n-l;
end