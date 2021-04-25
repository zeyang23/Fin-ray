% 有3根柔性支撑的finray，作用均匀压力

% 迭代逼近


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

Lcon1=3/4*sqrt((xA-xB)^2+(yA-yB)^2);
lambdaA1=1/4*LA;
lambdaB1=1/4*LB;

Lcon2=1/2*sqrt((xA-xB)^2+(yA-yB)^2);
lambdaA2=1/2*LA;
lambdaB2=1/2*LB;

Lcon3=1/4*sqrt((xA-xB)^2+(yA-yB)^2);
lambdaA3=3/4*LA;
lambdaB3=3/4*LB;


lambda1=0.4*LA;
lambda2=0.6*LA;


press_s=10;
press_e=30;

press_series=linspace(press_s,press_e,20);


problem_info=struct();
problem_info.Lcon1=Lcon1;
problem_info.Lcon2=Lcon2;
problem_info.Lcon3=Lcon3;

problem_info.lambdaA1=lambdaA1;
problem_info.lambdaA2=lambdaA2;
problem_info.lambdaA3=lambdaA3;

problem_info.lambdaB1=lambdaB1;
problem_info.lambdaB2=lambdaB2;
problem_info.lambdaB3=lambdaB3;

problem_info.lambda1=lambda1;
problem_info.lambda2=lambda2;


x0=zeros(27,1);
x_series=[];

for i=1:length(press_series)
    
    problem_info.press=press_series(i);


    f=@(x) check_balance(x,problem_info,psi,pA,pB,LA,LB,EA,EB,EC,IA,IB,IC);

    options_A = optimoptions('fsolve','Algorithm','levenberg-marquardt');
    options_B = optimoptions('fsolve','Algorithm','trust-region');
    options_C = optimoptions('fsolve','Algorithm','trust-region-dogleg');
    [x_series(:,i),fval_cosserat,exitflag_cosserat,output_cosserat]=fsolve(f,x0,options_B);
    
    if exitflag_cosserat<=0
        error('failed')
    end
    
    x0=x_series(:,i);
end



%% 验证求解结果
x_cosserat=x_series(:,end);

press=press_series(end);
problem_info.press=press;

na_0=x_cosserat(1:2);
ma_0=x_cosserat(3);

nb_0=x_cosserat(4:5);
mb_0=x_cosserat(6);

pe=x_cosserat(7:9);

f_p1=x_cosserat(10:11);
l_p1=x_cosserat(12);

f_q1=x_cosserat(13:14);
l_q1=x_cosserat(15);

f_p2=x_cosserat(16:17);
l_p2=x_cosserat(18);

f_q2=x_cosserat(19:20);
l_q2=x_cosserat(21);

f_p3=x_cosserat(22:23);
l_p3=x_cosserat(24);

f_q3=x_cosserat(25:26);
l_q3=x_cosserat(27);
    
    
options_a=odeset('MaxStep',1e-2);

span_a1 = [0 lambdaA1];
ya1_0 = [pA;na_0;ma_0];

[sA_1,YA_1] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a1, ya1_0,options_a);
ya_1=transpose(YA_1(end,:));


span_a2 = [lambdaA1,lambdaA2];
ya2_0=ya_1+[0;0;0;-f_p1;-l_p1];
[sA_2,YA_2] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a2, ya2_0,options_a);
ya_2=transpose(YA_2(end,:));

span_a3 = [lambdaA2,lambdaA3];
ya3_0=ya_2+[0;0;0;-f_p2;-l_p2];
[sA_3,YA_3] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a3, ya3_0,options_a);
ya_3=transpose(YA_3(end,:));

span_a4 = [lambdaA3,LA];
ya4_0=ya_3+[0;0;0;-f_p3;-l_p3];
[sA_4,YA_4] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a4, ya4_0,options_a);
ya_4=transpose(YA_4(end,:));


options_b=odeset('MaxStep',1e-2);

span_b1 = [0 lambdaB1];
yb1_0 = [pB;nb_0;mb_0];

[sB_1,YB_1] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b1, yb1_0,options_b);
yb_1=transpose(YB_1(end,:));


span_b2 = [lambdaB1,lambdaB2];
yb2_0=yb_1+[0;0;0;-f_q1;-l_q1];
[sB_2,YB_2] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b2, yb2_0,options_b);
yb_2=transpose(YB_2(end,:));

span_b3 = [lambdaB2,lambdaB3];
yb3_0=yb_2+[0;0;0;-f_q2;-l_q2];
[sB_3,YB_3] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b3, yb3_0,options_b);
yb_3=transpose(YB_3(end,:));

span_b4 = [lambdaB3,LB];
yb4_0=yb_3+[0;0;0;-f_q3;-l_q3];
[sB_4,YB_4] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b4, yb4_0,options_b);
yb_4=transpose(YB_4(end,:));


alpha=pi/2-psi/2;
beta=pi-alpha;

options_c1=odeset('MaxStep',1e-2);
span_c1=[0,Lcon1];
yc1_0=[ya_1(1:2);ya_1(3)-alpha;f_p1;l_p1];
[~,YC_1] = ode45(@(s,y) get_ydot_B(s,y,Lcon1,EC,IC), span_c1, yc1_0,options_c1);
ye_c1=transpose(YC_1(end,:));

options_c2=odeset('MaxStep',1e-2);
span_c2=[0,Lcon2];
yc2_0=[ya_2(1:2);ya_2(3)-alpha;f_p2;l_p2];
[~,YC_2] = ode45(@(s,y) get_ydot_B(s,y,Lcon2,EC,IC), span_c2, yc2_0,options_c2);
ye_c2=transpose(YC_2(end,:));

options_c3=odeset('MaxStep',1e-2);
span_c3=[0,Lcon3];
yc3_0=[ya_3(1:2);ya_3(3)-alpha;f_p3;l_p3];
[~,YC_3] = ode45(@(s,y) get_ydot_B(s,y,Lcon3,EC,IC), span_c3, yc3_0,options_c3);
ye_c3=transpose(YC_3(end,:));

hold on

plot(YC_1(:,1),YC_1(:,2))
plot(YC_2(:,1),YC_2(:,2))
plot(YC_3(:,1),YC_3(:,2))
    

YA_total=[YA_1;YA_2(2:end,:);YA_3(2:end,:);YA_4(2:end,:)];
YB_total=[YB_1;YB_2(2:end,:);YB_3(2:end,:);YB_4(2:end,:)];

sA_total=[sA_1;sA_2(2:end,:);sA_3(2:end,:);sA_4(2:end,:)];
sB_total=[sB_1;sB_2(2:end,:);sB_3(2:end,:);sB_4(2:end,:)];


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
function res=check_balance(x,problem_info,psi,pA,pB,LA,LB,EA,EB,EC,IA,IB,IC)
    Lcon1=problem_info.Lcon1;
    Lcon2=problem_info.Lcon2;
    Lcon3=problem_info.Lcon3;

    lambdaA1=problem_info.lambdaA1;
    lambdaA2=problem_info.lambdaA2;
    lambdaA3=problem_info.lambdaA3;

    lambdaB1=problem_info.lambdaB1;
    lambdaB2=problem_info.lambdaB2;
    lambdaB3=problem_info.lambdaB3;

    lambda1=problem_info.lambda1;
    lambda2=problem_info.lambda2;
    press=problem_info.press;
    
    
    na_0=x(1:2);
    ma_0=x(3);
    
    nb_0=x(4:5);
    mb_0=x(6);
    
    pe=x(7:9);
    
    f_p1=x(10:11);
    l_p1=x(12);
    
    f_q1=x(13:14);
    l_q1=x(15);
    
    f_p2=x(16:17);
    l_p2=x(18);
    
    f_q2=x(19:20);
    l_q2=x(21);
    
    f_p3=x(22:23);
    l_p3=x(24);
    
    f_q3=x(25:26);
    l_q3=x(27);
    
    
    options_a=odeset('MaxStep',1e-2);
    
    span_a1 = [0 lambdaA1];
    ya1_0 = [pA;na_0;ma_0];
    
    [~,YA_1] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a1, ya1_0,options_a);
    ya_1=transpose(YA_1(end,:));
    
    
    span_a2 = [lambdaA1,lambdaA2];
    ya2_0=ya_1+[0;0;0;-f_p1;-l_p1];
    [~,YA_2] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a2, ya2_0,options_a);
    ya_2=transpose(YA_2(end,:));
    
    span_a3 = [lambdaA2,lambdaA3];
    ya3_0=ya_2+[0;0;0;-f_p2;-l_p2];
    [~,YA_3] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a3, ya3_0,options_a);
    ya_3=transpose(YA_3(end,:));
    
    span_a4 = [lambdaA3,LA];
    ya4_0=ya_3+[0;0;0;-f_p3;-l_p3];
    [~,YA_4] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a4, ya4_0,options_a);
    ya_4=transpose(YA_4(end,:));
    
    
    options_b=odeset('MaxStep',1e-2);
    
    span_b1 = [0 lambdaB1];
    yb1_0 = [pB;nb_0;mb_0];
    
    [~,YB_1] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b1, yb1_0,options_b);
    yb_1=transpose(YB_1(end,:));
    
    
    span_b2 = [lambdaB1,lambdaB2];
    yb2_0=yb_1+[0;0;0;-f_q1;-l_q1];
    [~,YB_2] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b2, yb2_0,options_b);
    yb_2=transpose(YB_2(end,:));
    
    span_b3 = [lambdaB2,lambdaB3];
    yb3_0=yb_2+[0;0;0;-f_q2;-l_q2];
    [~,YB_3] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b3, yb3_0,options_b);
    yb_3=transpose(YB_3(end,:));
    
    span_b4 = [lambdaB3,LB];
    yb4_0=yb_3+[0;0;0;-f_q3;-l_q3];
    [~,YB_4] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b4, yb4_0,options_b);
    yb_4=transpose(YB_4(end,:));
    
    
    alpha=pi/2-psi/2;
    beta=pi-alpha;
    
    options_c1=odeset('MaxStep',1e-2);
    span_c1=[0,Lcon1];
    yc1_0=[ya_1(1:2);ya_1(3)-alpha;f_p1;l_p1];
    [~,YC_1] = ode45(@(s,y) get_ydot_B(s,y,Lcon1,EC,IC), span_c1, yc1_0,options_c1);
    ye_c1=transpose(YC_1(end,:));
    
    options_c2=odeset('MaxStep',1e-2);
    span_c2=[0,Lcon2];
    yc2_0=[ya_2(1:2);ya_2(3)-alpha;f_p2;l_p2];
    [~,YC_2] = ode45(@(s,y) get_ydot_B(s,y,Lcon2,EC,IC), span_c2, yc2_0,options_c2);
    ye_c2=transpose(YC_2(end,:));
    
    options_c3=odeset('MaxStep',1e-2);
    span_c3=[0,Lcon3];
    yc3_0=[ya_3(1:2);ya_3(3)-alpha;f_p3;l_p3];
    [~,YC_3] = ode45(@(s,y) get_ydot_B(s,y,Lcon3,EC,IC), span_c3, yc3_0,options_c3);
    ye_c3=transpose(YC_3(end,:));
    
    
    res=zeros(27,1);
    res(1:2)=ya_4(1:2)-pe(1:2);
    res(3)=ya_4(3)-pe(3);
    
    res(4:5)=yb_4(1:2)-pe(1:2);
    res(6)=yb_4(3)-pe(3)-psi;
    
    res(7:8)=ya_4(4:5)+yb_4(4:5);
    res(9)=transpose(ya_4(1:2))*[0 1;-1 0]*ya_4(4:5)+ya_4(6)...
          +transpose(yb_4(1:2))*[0 1;-1 0]*yb_4(4:5)+yb_4(6);
    
    res(10:11)=ye_c1(4:5)+f_q1;
    res(12)=ye_c1(6)+l_q1;
    
    res(13:14)=ye_c1(1:2)-yb_1(1:2);
    res(15)=ye_c1(3)+beta-yb_1(3);
    
    
    res(16:17)=ye_c2(4:5)+f_q2;
    res(18)=ye_c2(6)+l_q2;
    
    res(19:20)=ye_c2(1:2)-yb_2(1:2);
    res(21)=ye_c2(3)+beta-yb_2(3);
    
    
    res(22:23)=ye_c3(4:5)+f_q3;
    res(24)=ye_c3(6)+l_q3;
    
    res(25:26)=ye_c3(1:2)-yb_3(1:2);
    res(27)=ye_c3(3)+beta-yb_3(3);

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