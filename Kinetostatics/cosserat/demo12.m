% 尝试根据曲率来反解三角函数载荷的大小，载荷位置已知

%% 正向计算，生成数据

clear
clc

wid=14e-3;
thi=0.3e-3;
E=197*1e9;

Iz=1/12*thi.^3*wid;


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


Lcon1=3/4*sqrt((xA-xB)^2+(yA-yB)^2);
lambdaA1=1/4*LA;
lambdaB1=1/4*LB;

Lcon2=1/2*sqrt((xA-xB)^2+(yA-yB)^2);
lambdaA2=1/2*LA;
lambdaB2=1/2*LB;

Lcon3=1/4*sqrt((xA-xB)^2+(yA-yB)^2);
lambdaA3=3/4*LA;
lambdaB3=3/4*LB;


lambda1=0.2*LA;
lambda2=0.8*LA;

press=0.5;


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
problem_info.press=press;

f=@(x) check_balance(x,problem_info,psi,pA,pB,LA,LB,EA,EB,IA,IB);


x0=zeros(15,1);

options_A = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
options_B = optimoptions('fsolve','Display','off','Algorithm','trust-region');
options_C = optimoptions('fsolve','Display','off','Algorithm','trust-region-dogleg');
[x_cosserat,fval_cosserat,exitflag_cosserat,output_cosserat]=fsolve(f,x0,options_C);
% [x_cosserat,resnorm,residual,exitflag_cosserat,output_cosserat] = lsqnonlin(f,x0);


% 画出正向求解结果
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

[sA_1,YA_1] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a1, ya1_0,options_a);
ya_1=transpose(YA_1(end,:));


FA_1=[-fcon1*cos(gamma1);-fcon1*sin(gamma1)];
span_a2 = [lambdaA1,lambdaA2];
ya2_0=ya_1+[0;0;0;-FA_1;0];
[sA_2,YA_2] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a2, ya2_0,options_a);
ya_2=transpose(YA_2(end,:));

FA_2=[-fcon2*cos(gamma2);-fcon2*sin(gamma2)];
span_a3 = [lambdaA2,lambdaA3];
ya3_0=ya_2+[0;0;0;-FA_2;0];
[sA_3,YA_3] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a3, ya3_0,options_a);
ya_3=transpose(YA_3(end,:));

FA_3=[-fcon3*cos(gamma3);-fcon3*sin(gamma3)];
span_a4 = [lambdaA3,LA];
ya4_0=ya_3+[0;0;0;-FA_3;0];
[sA_4,YA_4] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a4, ya4_0,options_a);
ya_4=transpose(YA_4(end,:));


options_b=odeset('MaxStep',1e-2);

span_b1 = [0 lambdaB1];
yb1_0 = [pB;nb_0;mb_0];

[sB_1,YB_1] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b1, yb1_0,options_b);
yb_1=transpose(YB_1(end,:));


FB_1=[fcon1*cos(gamma1);fcon1*sin(gamma1)];
span_b2 = [lambdaB1,lambdaB2];
yb2_0=yb_1+[0;0;0;-FB_1;0];
[sB_2,YB_2] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b2, yb2_0,options_b);
yb_2=transpose(YB_2(end,:));

FB_2=[fcon2*cos(gamma2);fcon2*sin(gamma2)];
span_b3 = [lambdaB2,lambdaB3];
yb3_0=yb_2+[0;0;0;-FB_2;0];
[sB_3,YB_3] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b3, yb3_0,options_b);
yb_3=transpose(YB_3(end,:));

FB_3=[fcon3*cos(gamma3);fcon3*sin(gamma3)];
span_b4 = [lambdaB3,LB];
yb4_0=yb_3+[0;0;0;-FB_3;0];
[sB_4,YB_4] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b4, yb4_0,options_b);
yb_4=transpose(YB_4(end,:));



YA_total=[YA_1;YA_2(2:end,:);YA_3(2:end,:);YA_4(2:end,:)];
YB_total=[YB_1;YB_2(2:end,:);YB_3(2:end,:);YB_4(2:end,:)];

sA_total=[sA_1;sA_2(2:end,:);sA_3(2:end,:);sA_4(2:end,:)];
sB_total=[sB_1;sB_2(2:end,:);sB_3(2:end,:);sB_4(2:end,:)];

hold on
plot(YA_total(:,1),YA_total(:,2))
plot(YB_total(:,1),YB_total(:,2))
axis equal

plot([ya_1(1) yb_1(1)],[ya_1(2) yb_1(2)])
plot([ya_2(1) yb_2(1)],[ya_2(2) yb_2(2)])
plot([ya_3(1) yb_3(1)],[ya_3(2) yb_3(2)])

press_start_x=interp1(sA_total,YA_total(:,1),lambda1);
press_start_y=interp1(sA_total,YA_total(:,2),lambda1);
press_end_x=interp1(sA_total,YA_total(:,1),lambda2);
press_end_y=interp1(sA_total,YA_total(:,2),lambda2);

plot(press_start_x,press_start_y,'o')
plot(press_end_x,press_end_y,'o')


lambda_K=0.5*LB;

M1=interp1(sB_total,YB_total(:,6),lambda_K);
K1_true=M1/EB/IB;

% 添加误差
noise=0.2*(-1 + 2*rand);
K1=K1_true+noise;

fprintf('K1_true: %f\n',K1_true);
fprintf('K1: %f\n',K1);

fprintf('press_true: %f\n',press);



%% 根据曲率反解柔性板形状
problem_info_shape=struct();
problem_info_shape.Lcon1=Lcon1;
problem_info_shape.Lcon2=Lcon2;
problem_info_shape.Lcon3=Lcon3;

problem_info_shape.lambdaA1=lambdaA1;
problem_info_shape.lambdaA2=lambdaA2;
problem_info_shape.lambdaA3=lambdaA3;

problem_info_shape.lambdaB1=lambdaB1;
problem_info_shape.lambdaB2=lambdaB2;
problem_info_shape.lambdaB3=lambdaB3;

problem_info_shape.lambda1=lambda1;
problem_info_shape.lambda2=lambda2;

problem_info_shape.lambda_K=lambda_K;
problem_info_shape.K1=K1;


g=@(x) check_shape(x,problem_info_shape,psi,pA,pB,LA,LB,EA,EB,IA,IB);


x0_shape=zeros(16,1);

options_A = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
options_B = optimoptions('fsolve','Display','off','Algorithm','trust-region');
options_C = optimoptions('fsolve','Display','off','Algorithm','trust-region-dogleg');
[x_shape,fval_shape,exitflag_shape,output_shape]=fsolve(g,x0_shape,options_C);
% [x_shape,resnorm_shape,residual_shape,exitflag_shape,output_shape] = lsqnonlin(g,x0_shape);

save('x_shape.mat','x_shape');

%% 画出形状反解结果
clear

load('x_shape.mat');

fprintf('press_true: %f\n',x_shape(16));

wid=14e-3;
thi=0.3e-3;
E=197*1e9;

Iz=1/12*thi.^3*wid;


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


Lcon1=3/4*sqrt((xA-xB)^2+(yA-yB)^2);
lambdaA1=1/4*LA;
lambdaB1=1/4*LB;

Lcon2=1/2*sqrt((xA-xB)^2+(yA-yB)^2);
lambdaA2=1/2*LA;
lambdaB2=1/2*LB;

Lcon3=1/4*sqrt((xA-xB)^2+(yA-yB)^2);
lambdaA3=3/4*LA;
lambdaB3=3/4*LB;


lambda1=0.2*LA;
lambda2=0.8*LA;



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
problem_info.press=x_shape(16);

f=@(x) check_balance(x,problem_info,psi,pA,pB,LA,LB,EA,EB,IA,IB);


x0=zeros(15,1);

options_A = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
options_B = optimoptions('fsolve','Display','off','Algorithm','trust-region');
options_C = optimoptions('fsolve','Display','off','Algorithm','trust-region-dogleg');
[x_cosserat,fval_cosserat,exitflag_cosserat,output_cosserat]=fsolve(f,x0,options_C);
% [x_cosserat,resnorm,residual,exitflag_cosserat,output_cosserat] = lsqnonlin(f,x0);


% 画出求解结果
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

press=problem_info.press;


options_a=odeset('MaxStep',1e-2);

span_a1 = [0 lambdaA1];
ya1_0 = [pA;na_0;ma_0];

[sA_1,YA_1] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a1, ya1_0,options_a);
ya_1=transpose(YA_1(end,:));


FA_1=[-fcon1*cos(gamma1);-fcon1*sin(gamma1)];
span_a2 = [lambdaA1,lambdaA2];
ya2_0=ya_1+[0;0;0;-FA_1;0];
[sA_2,YA_2] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a2, ya2_0,options_a);
ya_2=transpose(YA_2(end,:));

FA_2=[-fcon2*cos(gamma2);-fcon2*sin(gamma2)];
span_a3 = [lambdaA2,lambdaA3];
ya3_0=ya_2+[0;0;0;-FA_2;0];
[sA_3,YA_3] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a3, ya3_0,options_a);
ya_3=transpose(YA_3(end,:));

FA_3=[-fcon3*cos(gamma3);-fcon3*sin(gamma3)];
span_a4 = [lambdaA3,LA];
ya4_0=ya_3+[0;0;0;-FA_3;0];
[sA_4,YA_4] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a4, ya4_0,options_a);
ya_4=transpose(YA_4(end,:));


options_b=odeset('MaxStep',1e-2);

span_b1 = [0 lambdaB1];
yb1_0 = [pB;nb_0;mb_0];

[sB_1,YB_1] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b1, yb1_0,options_b);
yb_1=transpose(YB_1(end,:));


FB_1=[fcon1*cos(gamma1);fcon1*sin(gamma1)];
span_b2 = [lambdaB1,lambdaB2];
yb2_0=yb_1+[0;0;0;-FB_1;0];
[sB_2,YB_2] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b2, yb2_0,options_b);
yb_2=transpose(YB_2(end,:));

FB_2=[fcon2*cos(gamma2);fcon2*sin(gamma2)];
span_b3 = [lambdaB2,lambdaB3];
yb3_0=yb_2+[0;0;0;-FB_2;0];
[sB_3,YB_3] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b3, yb3_0,options_b);
yb_3=transpose(YB_3(end,:));

FB_3=[fcon3*cos(gamma3);fcon3*sin(gamma3)];
span_b4 = [lambdaB3,LB];
yb4_0=yb_3+[0;0;0;-FB_3;0];
[sB_4,YB_4] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b4, yb4_0,options_b);
yb_4=transpose(YB_4(end,:));



YA_total=[YA_1;YA_2(2:end,:);YA_3(2:end,:);YA_4(2:end,:)];
YB_total=[YB_1;YB_2(2:end,:);YB_3(2:end,:);YB_4(2:end,:)];

sA_total=[sA_1;sA_2(2:end,:);sA_3(2:end,:);sA_4(2:end,:)];
sB_total=[sB_1;sB_2(2:end,:);sB_3(2:end,:);sB_4(2:end,:)];

hold on
plot(YA_total(:,1),YA_total(:,2))
plot(YB_total(:,1),YB_total(:,2))
axis equal

plot([ya_1(1) yb_1(1)],[ya_1(2) yb_1(2)])
plot([ya_2(1) yb_2(1)],[ya_2(2) yb_2(2)])
plot([ya_3(1) yb_3(1)],[ya_3(2) yb_3(2)])

press_start_x=interp1(sA_total,YA_total(:,1),lambda1);
press_start_y=interp1(sA_total,YA_total(:,2),lambda1);
press_end_x=interp1(sA_total,YA_total(:,1),lambda2);
press_end_y=interp1(sA_total,YA_total(:,2),lambda2);

plot(press_start_x,press_start_y,'o')
plot(press_end_x,press_end_y,'o')

%%
function res=check_shape(x,problem_info,psi,pA,pB,LA,LB,EA,EB,IA,IB)
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
    
    lambda_K=problem_info.lambda_K;
    K1=problem_info.K1;
    
    
    press=x(16);
    
    
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
    
    [sA_1,YA_1] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a1, ya1_0,options_a);
    ya_1=transpose(YA_1(end,:));
    
    
    FA_1=[-fcon1*cos(gamma1);-fcon1*sin(gamma1)];
    span_a2 = [lambdaA1,lambdaA2];
    ya2_0=ya_1+[0;0;0;-FA_1;0];
    [sA_2,YA_2] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a2, ya2_0,options_a);
    ya_2=transpose(YA_2(end,:));
    
    FA_2=[-fcon2*cos(gamma2);-fcon2*sin(gamma2)];
    span_a3 = [lambdaA2,lambdaA3];
    ya3_0=ya_2+[0;0;0;-FA_2;0];
    [sA_3,YA_3] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a3, ya3_0,options_a);
    ya_3=transpose(YA_3(end,:));
    
    FA_3=[-fcon3*cos(gamma3);-fcon3*sin(gamma3)];
    span_a4 = [lambdaA3,LA];
    ya4_0=ya_3+[0;0;0;-FA_3;0];
    [sA_4,YA_4] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a4, ya4_0,options_a);
    ya_4=transpose(YA_4(end,:));
    
    
    options_b=odeset('MaxStep',1e-2);
    
    span_b1 = [0 lambdaB1];
    yb1_0 = [pB;nb_0;mb_0];
    
    [sB_1,YB_1] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b1, yb1_0,options_b);
    yb_1=transpose(YB_1(end,:));
    
    
    FB_1=[fcon1*cos(gamma1);fcon1*sin(gamma1)];
    span_b2 = [lambdaB1,lambdaB2];
    yb2_0=yb_1+[0;0;0;-FB_1;0];
    [sB_2,YB_2] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b2, yb2_0,options_b);
    yb_2=transpose(YB_2(end,:));
    
    FB_2=[fcon2*cos(gamma2);fcon2*sin(gamma2)];
    span_b3 = [lambdaB2,lambdaB3];
    yb3_0=yb_2+[0;0;0;-FB_2;0];
    [sB_3,YB_3] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b3, yb3_0,options_b);
    yb_3=transpose(YB_3(end,:));
    
    FB_3=[fcon3*cos(gamma3);fcon3*sin(gamma3)];
    span_b4 = [lambdaB3,LB];
    yb4_0=yb_3+[0;0;0;-FB_3;0];
    [sB_4,YB_4] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b4, yb4_0,options_b);
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
    
    
    YA_total=[YA_1;YA_2(2:end,:);YA_3(2:end,:);YA_4(2:end,:)];
    YB_total=[YB_1;YB_2(2:end,:);YB_3(2:end,:);YB_4(2:end,:)];

    sA_total=[sA_1;sA_2(2:end,:);sA_3(2:end,:);sA_4(2:end,:)];
    sB_total=[sB_1;sB_2(2:end,:);sB_3(2:end,:);sB_4(2:end,:)];
    
    M1_now=interp1(sB_total,YB_total(:,6),lambda_K);
    K1_now=M1_now/EB/IB;
    
    res(16)=K1_now-K1;
    
    
end

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

    lambda1=problem_info.lambda1;
    lambda2=problem_info.lambda2;
    press=problem_info.press;
    
    
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
    
    [~,YA_1] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a1, ya1_0,options_a);
    ya_1=transpose(YA_1(end,:));
    
    
    FA_1=[-fcon1*cos(gamma1);-fcon1*sin(gamma1)];
    span_a2 = [lambdaA1,lambdaA2];
    ya2_0=ya_1+[0;0;0;-FA_1;0];
    [~,YA_2] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a2, ya2_0,options_a);
    ya_2=transpose(YA_2(end,:));
    
    FA_2=[-fcon2*cos(gamma2);-fcon2*sin(gamma2)];
    span_a3 = [lambdaA2,lambdaA3];
    ya3_0=ya_2+[0;0;0;-FA_2;0];
    [~,YA_3] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a3, ya3_0,options_a);
    ya_3=transpose(YA_3(end,:));
    
    FA_3=[-fcon3*cos(gamma3);-fcon3*sin(gamma3)];
    span_a4 = [lambdaA3,LA];
    ya4_0=ya_3+[0;0;0;-FA_3;0];
    [~,YA_4] = ode45(@(s,y) get_ydot_A(s,y,lambda1,lambda2,press,LA,EA,IA), span_a4, ya4_0,options_a);
    ya_4=transpose(YA_4(end,:));
    
    
    options_b=odeset('MaxStep',1e-2);
    
    span_b1 = [0 lambdaB1];
    yb1_0 = [pB;nb_0;mb_0];
    
    [~,YB_1] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b1, yb1_0,options_b);
    yb_1=transpose(YB_1(end,:));
    
    
    FB_1=[fcon1*cos(gamma1);fcon1*sin(gamma1)];
    span_b2 = [lambdaB1,lambdaB2];
    yb2_0=yb_1+[0;0;0;-FB_1;0];
    [~,YB_2] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b2, yb2_0,options_b);
    yb_2=transpose(YB_2(end,:));
    
    FB_2=[fcon2*cos(gamma2);fcon2*sin(gamma2)];
    span_b3 = [lambdaB2,lambdaB3];
    yb3_0=yb_2+[0;0;0;-FB_2;0];
    [~,YB_3] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b3, yb3_0,options_b);
    yb_3=transpose(YB_3(end,:));
    
    FB_3=[fcon3*cos(gamma3);fcon3*sin(gamma3)];
    span_b4 = [lambdaB3,LB];
    yb4_0=yb_3+[0;0;0;-FB_3;0];
    [~,YB_4] = ode45(@(s,y) get_ydot_B(s,y,LB,EB,IB), span_b4, yb4_0,options_b);
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
        f=press*sin((s-lambda1)/(lambda2-lambda1)*pi);
        fs=[f*sin(theta);-f*cos(theta)];
    else
        fs=[0;0];
    end
end