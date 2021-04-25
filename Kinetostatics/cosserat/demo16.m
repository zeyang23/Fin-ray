% 有外力、柔性支撑的finray
% 逐渐迭代求解，为了解决初值敏感性的问题



%% Cosserat Model shooting method
clear
clc

wid=14e-3;
thi=0.3e-3;
E=197*1e9;

Iz=1/12*thi.^3*wid;


L0=0.25;

LA=L0;
LB=L0;

EA=E;
EB=E;

EC=0.5e9;

IA=Iz;
IB=Iz;
IC=Iz;


half_psi_degree=10;
psi_degree=20;


alpha_degree=90-half_psi_degree;
beta_degree=180-alpha_degree;

psi=psi_degree/180*pi;
alpha=alpha_degree/180*pi;
beta=beta_degree/180*pi;


xA=0;
yA=0;

xB=2*L0*sin(psi/2);
yB=0;


pA=[xA;yA;alpha];
pB=[xB;yB;beta];

Lcon=1/2*sqrt((xA-xB)^2+(yA-yB)^2);
lambdaA=1/2*LA;
lambdaB=1/2*LB;

% 定义均匀载荷
lambda1=0.4*LA;
lambda2=0.6*LA;

press_s=25;
press_e=55;

press_series=linspace(press_s,press_e,10);


x0=zeros(15,1);

% x0(7:9)=[L0*sin(psi/2);L0*cos(psi/2);alpha];


options_A = optimoptions('fsolve','Algorithm','levenberg-marquardt');
options_B = optimoptions('fsolve','Algorithm','trust-region');
options_C = optimoptions('fsolve','Algorithm','trust-region-dogleg');

x_series=[];

for i=1:length(press_series)
    press=press_series(i);
    
    f=@(x) check_balance(x,lambda1,lambda2,press,Lcon,lambdaA,lambdaB,psi,pA,pB,LA,LB,EA,EB,EC,IA,IB,IC);
    
    [x_series(:,i),fval_cosserat,exitflag_cosserat,output_cosserat]=fsolve(f,x0,options_B);
    
    if exitflag_cosserat<=0
        error('failed')
    end
    
    x0=x_series(:,i);
end

%% 验证求解结果
x_cosserat=x_series(:,end);

na_0=x_cosserat(1:2);
ma_0=x_cosserat(3);

nb_0=x_cosserat(4:5);
mb_0=x_cosserat(6);

pe=x_cosserat(7:9);

f_p=x_cosserat(10:11);
l_p=x_cosserat(12);

f_q=x_cosserat(13:14);
l_q=x_cosserat(15);
    
    
span_a1 = [0 lambdaA];
ya1_0 = [pA;na_0;ma_0];
options_a=odeset('MaxStep',1e-2);
[sA_in,YA_in] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a1, ya1_0,options_a);
yin_a=transpose(YA_in(end,:));


span_a2 = [lambdaA,LA];
ya2_0=yin_a+[0;0;0;-f_p;-l_p];
[sA,YA] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a2, ya2_0,options_a);
ye_a=transpose(YA(end,:));


span_b1 = [0 lambdaB];
yb1_0 = [pB;nb_0;mb_0];
options_b=odeset('MaxStep',1e-2);
[sB_in,YB_in] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b1, yb1_0,options_b);
yin_b=transpose(YB_in(end,:));


span_b2 = [lambdaB,LB];
yb2_0=yin_b+[0;0;0;-f_q;-l_q];
[sB,YB] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b2, yb2_0,options_b);
ye_b=transpose(YB(end,:));


options_c=odeset('MaxStep',1e-2);
span_c=[0,Lcon];
alpha=pi/2-psi/2;
beta=pi-alpha;
yc_0=[yin_a(1:2);yin_a(3)-alpha;f_p;l_p];
[~,YC] = ode45(@(s,y) get_ydot_B(s,y,Lcon,EC,IC), span_c, yc_0,options_c);


plot(YC(:,1),YC(:,2))
    

YA_total=[YA_in;YA(2:end,:)];
YB_total=[YB_in;YB(2:end,:)];

sA_total=[sA_in;sA(2:end,:)];
sB_total=[sB_in;sB(2:end,:)];

hold on
plot(YA_total(:,1),YA_total(:,2))
plot(YB_total(:,1),YB_total(:,2))
axis equal

% plot([yin_a(1) yin_b(1)],[yin_a(2) yin_b(2)])

press_start_x=interp1(sA_total,YA_total(:,1),lambda1);
press_start_y=interp1(sA_total,YA_total(:,2),lambda1);
press_end_x=interp1(sA_total,YA_total(:,1),lambda2);
press_end_y=interp1(sA_total,YA_total(:,2),lambda2);

plot(press_start_x,press_start_y,'o')
plot(press_end_x,press_end_y,'o')


%%
function res=check_balance(x,lambda1,lambda2,press,Lcon,lambdaA,lambdaB,psi,pA,pB,LA,LB,EA,EB,EC,IA,IB,IC)
    
    na_0=x(1:2);
    ma_0=x(3);
    
    nb_0=x(4:5);
    mb_0=x(6);
    
    pe=x(7:9);
    
    f_p=x(10:11);
    l_p=x(12);
    
    f_q=x(13:14);
    l_q=x(15);
    
    
    span_a1 = [0 lambdaA];
    ya1_0 = [pA;na_0;ma_0];
    options_a=odeset('MaxStep',1e-2);
    [~,YA_in] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a1, ya1_0,options_a);
    yin_a=transpose(YA_in(end,:));
    
    
    span_a2 = [lambdaA,LA];
    ya2_0=yin_a+[0;0;0;-f_p;-l_p];
    [~,YA] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a2, ya2_0,options_a);
    ye_a=transpose(YA(end,:));
    
    
    span_b1 = [0 lambdaB];
    yb1_0 = [pB;nb_0;mb_0];
    options_b=odeset('MaxStep',1e-2);
    [~,YB_in] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b1, yb1_0,options_b);
    yin_b=transpose(YB_in(end,:));
    
    
    span_b2 = [lambdaB,LB];
    yb2_0=yin_b+[0;0;0;-f_q;-l_q];
    [~,YB] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b2, yb2_0,options_b);
    ye_b=transpose(YB(end,:));
    
    
    options_c=odeset('MaxStep',1e-2);
    span_c=[0,Lcon];
    alpha=pi/2-psi/2;
    beta=pi-alpha;
    yc_0=[yin_a(1:2);yin_a(3)-alpha;f_p;l_p];
    [~,YC] = ode45(@(s,y) get_ydot_B(s,y,Lcon,EC,IC), span_c, yc_0,options_c);
    
    ye_c=transpose(YC(end,:));
    
    
    res=zeros(15,1);
    res(1:2)=ye_a(1:2)-pe(1:2);
    res(3)=ye_a(3)-pe(3);
    
    res(4:5)=ye_b(1:2)-pe(1:2);
    res(6)=ye_b(3)-pe(3)-psi;
    
    res(7:8)=ye_a(4:5)+ye_b(4:5);
    res(9)=transpose(ye_a(1:2))*[0 1;-1 0]*ye_a(4:5)+ye_a(6)...
          +transpose(ye_b(1:2))*[0 1;-1 0]*ye_b(4:5)+ye_b(6);
    
    res(10:11)=ye_c(4:5)+f_q;
    res(12)=ye_c(6)+l_q;
    
    res(13:14)=ye_c(1:2)-yin_b(1:2);
    res(15)=ye_c(3)+beta-yin_b(3);

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