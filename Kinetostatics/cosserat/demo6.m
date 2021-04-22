% 无刚性约束的Finray
% 作用均匀法向载荷

clear
clc


%% Cosserat Model shooting method



wid=14e-3;
thi=0.3e-3;
E=197*1e9;

Iz=1/12*thi.^3*wid;

pdes=[0;0;pi];

L0=1;

LA=L0;
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


% 定义均匀载荷
lambda1=0.4*LA;
lambda2=0.6*LA;

press=0.5;

f=@(x) check_balance(x,lambda1,lambda2,press,psi,pA,pB,LA,LB,EA,EB,IA,IB);


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

sol_YA = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a, ya_0,options_a);


span_b = [0 LB];
yb_0 = [pB;nb_0;mb_0];
options_b=odeset('MaxStep',1e-2);
sol_YB = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b, yb_0,options_b);


hold on
plot(sol_YA.y(1,:),sol_YA.y(2,:))
plot(sol_YB.y(1,:),sol_YB.y(2,:))
axis equal

press_start=deval(sol_YA,lambda1);
press_end=deval(sol_YA,lambda2);

plot(press_start(1),press_start(2),'o')
plot(press_end(1),press_end(2),'o')

%%
function res=check_balance(x,lambda1,lambda2,press,psi,pA,pB,LA,LB,EA,EB,IA,IB)
    
    na_0=x(1:2);
    ma_0=x(3);
    
    nb_0=x(4:5);
    mb_0=x(6);
    
    pe=x(7:9);
    
    
    span_a = [0 LA];
    ya_0 = [pA;na_0;ma_0];
    options_a=odeset('MaxStep',1e-2);
    
    [~,YA] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a, ya_0,options_a);
    ye_a=transpose(YA(end,:));
    
    
    span_b = [0 LB];
    yb_0 = [pB;nb_0;mb_0];
    options_b=odeset('MaxStep',1e-2);
    [~,YB] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b, yb_0,options_b);
    ye_b=transpose(YB(end,:));
    
    
    res=zeros(9,1);
    res(1:2)=ye_a(1:2)-pe(1:2);
    res(3)=ye_a(3)-pe(3);
    
    res(4:5)=ye_b(1:2)-pe(1:2);
    res(6)=ye_b(3)-pe(3)-psi;
    
    res(7:8)=ye_a(4:5)+ye_b(4:5);
    res(9)=transpose(ye_a(1:2))*[0 1;-1 0]*ye_a(4:5)+ye_a(6)...
          +transpose(ye_b(1:2))*[0 1;-1 0]*ye_b(4:5)+ye_b(6);

end


function ydot=get_ydot_A(s,y,lambda1,lambda2,press,L,E,I)
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
    
    f=f_press(s,y,lambda1,lambda2,press);
    l=0;

    ydot(1:2)=[cos(theta);sin(theta)];
    ydot(3)=1/(E*I)*m;
    ydot(4:5)=-f;
    ydot(6)=-[cos(theta) sin(theta)]*[0 1;-1 0]*n-l;
end

function ydot=get_ydot_B(s,y,L,E,I)
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

function fs=f_press(s,y,lambda1,lambda2,press)
    theta=y(3);
    if(s>=lambda1) && (s<=lambda2)
        fs=[press*sin(theta);-press*cos(theta)];
    else
        fs=[0;0];
    end
end