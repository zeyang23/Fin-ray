% 使用cosserat模型求解无外力、有2根刚性约束的Finray
% 使用分段的方法来处理

% 21-04-23
% 能够求解，但是鲁棒性存疑，左侧的长度如果太短，则无法求解


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


Lcon1=3/4*sqrt((xA-xB)^2+(yA-yB)^2);
lambdaA1=1/4*LA;
lambdaB1=1/4*LB;

Lcon2=1/2*sqrt((xA-xB)^2+(yA-yB)^2);
lambdaA2=1/2*LA;
lambdaB2=1/2*LB;

Lcon3=1/4*sqrt((xA-xB)^2+(yA-yB)^2);
lambdaA3=3/4*LA;
lambdaB3=3/4*LB;


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


f=@(x) check_balance(x,problem_info,psi,pA,pB,LA,LB,EA,EB,IA,IB);


x0=zeros(15,1);

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

fcon1=x_cosserat(10);
gamma1=x_cosserat(11);
    
fcon2=x_cosserat(12);
gamma2=x_cosserat(13);

fcon3=x_cosserat(14);
gamma3=x_cosserat(15);


options_a=odeset('MaxStep',1e-2);

span_a1 = [0 lambdaA1];
ya1_0 = [pA;na_0;ma_0];

[~,YA_1] = ode45(@(s,y) get_ydot(s,y,LA,EA,IA), span_a1, ya1_0,options_a);
ya_1=transpose(YA_1(end,:));


FA_1=[-fcon1*cos(gamma1);-fcon1*sin(gamma1)];
span_a2 = [lambdaA1,lambdaA2];
ya2_0=ya_1+[0;0;0;-FA_1;0];
[~,YA_2] = ode45(@(s,y) get_ydot(s,y,LA,EA,IA), span_a2, ya2_0,options_a);
ya_2=transpose(YA_2(end,:));

FA_2=[-fcon2*cos(gamma2);-fcon2*sin(gamma2)];
span_a3 = [lambdaA2,lambdaA3];
ya3_0=ya_2+[0;0;0;-FA_2;0];
[~,YA_3] = ode45(@(s,y) get_ydot(s,y,LA,EA,IA), span_a3, ya3_0,options_a);
ya_3=transpose(YA_3(end,:));

FA_3=[-fcon3*cos(gamma3);-fcon3*sin(gamma3)];
span_a4 = [lambdaA3,LA];
ya4_0=ya_3+[0;0;0;-FA_3;0];
[~,YA_4] = ode45(@(s,y) get_ydot(s,y,LA,EA,IA), span_a4, ya4_0,options_a);
ya_4=transpose(YA_4(end,:));


options_b=odeset('MaxStep',1e-2);

span_b1 = [0 lambdaB1];
yb1_0 = [pB;nb_0;mb_0];

[~,YB_1] = ode45(@(s,y) get_ydot(s,y,LB,EB,IB), span_b1, yb1_0,options_b);
yb_1=transpose(YB_1(end,:));


FB_1=[fcon1*cos(gamma1);fcon1*sin(gamma1)];
span_b2 = [lambdaB1,lambdaB2];
yb2_0=yb_1+[0;0;0;-FB_1;0];
[~,YB_2] = ode45(@(s,y) get_ydot(s,y,LB,EB,IB), span_b2, yb2_0,options_b);
yb_2=transpose(YB_2(end,:));

FB_2=[fcon2*cos(gamma2);fcon2*sin(gamma2)];
span_b3 = [lambdaB2,lambdaB3];
yb3_0=yb_2+[0;0;0;-FB_2;0];
[~,YB_3] = ode45(@(s,y) get_ydot(s,y,LB,EB,IB), span_b3, yb3_0,options_b);
yb_3=transpose(YB_3(end,:));

FB_3=[fcon3*cos(gamma3);fcon3*sin(gamma3)];
span_b4 = [lambdaB3,LB];
yb4_0=yb_3+[0;0;0;-FB_3;0];
[~,YB_4] = ode45(@(s,y) get_ydot(s,y,LB,EB,IB), span_b4, yb4_0,options_b);
yb_4=transpose(YB_4(end,:));



YA_total=[YA_1;YA_2(2:end,:);YA_3(2:end,:);YA_4(2:end,:)];
YB_total=[YB_1;YB_2(2:end,:);YB_3(2:end,:);YB_4(2:end,:)];

hold on
plot(YA_total(:,1),YA_total(:,2))
plot(YB_total(:,1),YB_total(:,2))
axis equal

plot([ya_1(1) yb_1(1)],[ya_1(2) yb_1(2)])
plot([ya_2(1) yb_2(1)],[ya_2(2) yb_2(2)])
plot([ya_3(1) yb_3(1)],[ya_3(2) yb_3(2)])



%%
function res=check_balance(x,problem_info,psi,pA,pB,LA,LB,EA,EB,IA,IB)
    Lcon1=problem_info.Lcon1;
    Lcon2=problem_info.Lcon2;
    Lcon3=problem_info.Lcon3;

    lambdaA1=problem_info.lambdaA1;
    lambdaA2=problem_info.lambdaA2;
    lambdaA3=problem_info.lambdaA3;

    lambdaB1=problem_info.lambdaB1;
    lambdaB2=problem_info.lambdaB2;
    lambdaB3=problem_info.lambdaB3;


    
    na_0=x(1:2);
    ma_0=x(3);
    
    nb_0=x(4:5);
    mb_0=x(6);
    
    pe=x(7:9);
    
    fcon1=x(10);
    gamma1=x(11);
    
    fcon2=x(12);
    gamma2=x(13);
    
    fcon3=x(14);
    gamma3=x(15);
    
    
    options_a=odeset('MaxStep',1e-2);
    
    span_a1 = [0 lambdaA1];
    ya1_0 = [pA;na_0;ma_0];
    
    [~,YA_1] = ode45(@(s,y) get_ydot(s,y,LA,EA,IA), span_a1, ya1_0,options_a);
    ya_1=transpose(YA_1(end,:));
    
    
    FA_1=[-fcon1*cos(gamma1);-fcon1*sin(gamma1)];
    span_a2 = [lambdaA1,lambdaA2];
    ya2_0=ya_1+[0;0;0;-FA_1;0];
    [~,YA_2] = ode45(@(s,y) get_ydot(s,y,LA,EA,IA), span_a2, ya2_0,options_a);
    ya_2=transpose(YA_2(end,:));
    
    FA_2=[-fcon2*cos(gamma2);-fcon2*sin(gamma2)];
    span_a3 = [lambdaA2,lambdaA3];
    ya3_0=ya_2+[0;0;0;-FA_2;0];
    [~,YA_3] = ode45(@(s,y) get_ydot(s,y,LA,EA,IA), span_a3, ya3_0,options_a);
    ya_3=transpose(YA_3(end,:));
    
    FA_3=[-fcon3*cos(gamma3);-fcon3*sin(gamma3)];
    span_a4 = [lambdaA3,LA];
    ya4_0=ya_3+[0;0;0;-FA_3;0];
    [~,YA_4] = ode45(@(s,y) get_ydot(s,y,LA,EA,IA), span_a4, ya4_0,options_a);
    ya_4=transpose(YA_4(end,:));
    
    
    options_b=odeset('MaxStep',1e-2);
    
    span_b1 = [0 lambdaB1];
    yb1_0 = [pB;nb_0;mb_0];
    
    [~,YB_1] = ode45(@(s,y) get_ydot(s,y,LB,EB,IB), span_b1, yb1_0,options_b);
    yb_1=transpose(YB_1(end,:));
    
    
    FB_1=[fcon1*cos(gamma1);fcon1*sin(gamma1)];
    span_b2 = [lambdaB1,lambdaB2];
    yb2_0=yb_1+[0;0;0;-FB_1;0];
    [~,YB_2] = ode45(@(s,y) get_ydot(s,y,LB,EB,IB), span_b2, yb2_0,options_b);
    yb_2=transpose(YB_2(end,:));
    
    FB_2=[fcon2*cos(gamma2);fcon2*sin(gamma2)];
    span_b3 = [lambdaB2,lambdaB3];
    yb3_0=yb_2+[0;0;0;-FB_2;0];
    [~,YB_3] = ode45(@(s,y) get_ydot(s,y,LB,EB,IB), span_b3, yb3_0,options_b);
    yb_3=transpose(YB_3(end,:));
    
    FB_3=[fcon3*cos(gamma3);fcon3*sin(gamma3)];
    span_b4 = [lambdaB3,LB];
    yb4_0=yb_3+[0;0;0;-FB_3;0];
    [~,YB_4] = ode45(@(s,y) get_ydot(s,y,LB,EB,IB), span_b4, yb4_0,options_b);
    yb_4=transpose(YB_4(end,:));
    
    
    res=zeros(15,1);
    res(1:2)=ya_4(1:2)-pe(1:2);
    res(3)=ya_4(3)-pe(3);
    
    res(4:5)=yb_4(1:2)-pe(1:2);
    res(6)=yb_4(3)-pe(3)-psi;
    
    res(7:8)=ya_4(4:5)+yb_4(4:5);
    res(9)=transpose(ya_4(1:2))*[0 1;-1 0]*ya_4(4:5)+ya_4(6)...
          +transpose(yb_4(1:2))*[0 1;-1 0]*yb_4(4:5)+yb_4(6);
      
    res(10:11)=ya_1(1:2)+[Lcon1*cos(gamma1);Lcon1*sin(gamma1)]-yb_1(1:2);
    res(12:13)=ya_2(1:2)+[Lcon2*cos(gamma2);Lcon2*sin(gamma2)]-yb_2(1:2);
    res(14:15)=ya_3(1:2)+[Lcon3*cos(gamma3);Lcon3*sin(gamma3)]-yb_3(1:2);

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