% 与椭圆多点接触，非均匀间隔


%%
clear
clc


a=30e-3;
b=40e-3;

center_x=0-a;
center_y=66.5e-3;

delta=25e-3;

f=@(tangent_var) myfunc(tangent_var,delta,center_x,center_y,a,b);



% 相切点的个数
N_tan=1;


lb=zeros(3*N_tan,1);
for i=1:N_tan
    lb(3*i-2)=-pi/2;
end

ub=zeros(3*N_tan,1);
for i=1:N_tan
    ub(3*i-2)=pi/2;
    ub(3*i-1)=1;
    ub(3*i)=10;
end

t_series=linspace(-pi/3,pi/3,N_tan);
ratio_series=linspace(0.2,0.8,N_tan);

X0=zeros(3*N_tan,1);
for i=1:N_tan
    X0(3*i-2)=t_series(i);
    X0(3*i-1)=ratio_series(i);
    X0(3*i)=1;
end


[X,resnorm,residual,exitflag_lsq,output_lsq] = lsqnonlin(f,X0,lb,ub);
% [X,fval_fsolve,exitflag_fsolve,output_fsolve]=fsolve(f,X0);


%% 验证结果
% 画出椭圆位置
funx=@(t) a*cos(t)+center_x;
funy=@(t) b*sin(t)+center_y;

tinterval=[-pi,pi];

fplot(funx,funy,tinterval)

axis equal;
hold on


xA=0-delta;
yA=0;

xB=73.24e-3-delta;
yB=0;


alpha_degree=90;
beta_degree=116;
psi_degree=beta_degree-alpha_degree;

psi=psi_degree/180*pi;
alpha=alpha_degree/180*pi;
beta=beta_degree/180*pi;


nA=50;
nB=50;


LA=(xB-xA)/tan(psi);
LB=(xB-xA)/sin(psi);

wid=28e-3;
thi=0.2e-3;
E=213e9;


constraint_ratio=[];

La_1=27e-3;
La_2=54e-3;
La_3=81e-3;
La_4=108e-3;

Lb_1=45.67e-3;
Lb_2=78.69e-3;
Lb_3=107.5e-3;
Lb_4=132.27e-3;

Lcon1=sqrt((LA-La_1)^2+(LB-Lb_1)^2-2*(LA-La_1)*(LB-Lb_1)*cos(psi));
Lcon2=sqrt((LA-La_2)^2+(LB-Lb_2)^2-2*(LA-La_2)*(LB-Lb_2)*cos(psi));
Lcon3=sqrt((LA-La_3)^2+(LB-Lb_3)^2-2*(LA-La_3)*(LB-Lb_3)*cos(psi));
Lcon4=sqrt((LA-La_4)^2+(LB-Lb_4)^2-2*(LA-La_4)*(LB-Lb_4)*cos(psi));

% constraint_ratio=[Lcon1,La_1/LA,Lb_1/LB;
%                   Lcon2,La_2/LA,Lb_2/LB;
%                   Lcon3,La_3/LA,Lb_3/LB;
%                   Lcon4,La_4/LA,Lb_4/LB];

    
A_force_ratio=[];

tangent_ratio=zeros(N_tan,1);
tangent_F=zeros(N_tan,1);
    
for i=1:N_tan
    tangent_ratio(i)=X(3*i-1);
    tangent_F(i)=X(3*i);
end

A_force_ratio=[tangent_F,tangent_ratio];


B_force_ratio=[];


% 根据上面设定的参数生成参数结构体
finray_info=struct();

finray_info.pA=[xA;yA;alpha];
finray_info.pB=[xB;yB;beta];

finray_info.LA=LA;
finray_info.LB=LB;

finray_info.nA=nA;
finray_info.nB=nB;

finray_info.psi=psi;

finray_info.wid_A=wid;
finray_info.wid_B=wid;

finray_info.thi_A=thi;
finray_info.thi_B=thi;

finray_info.E_A=LA/135e-3*E;
finray_info.E_B=LB/156e-3*E;

finray_info.constraint_ratio=constraint_ratio;

finray_info.A_force_ratio=A_force_ratio;

finray_info.B_force_ratio=B_force_ratio;


Finray1=finray_force(finray_info);




x0=zeros(nA+nB+3+2*Finray1.constraint_number,1);

options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

f=@(x) Finray1.cal_balance(x);
[x_solve,fval,exitflag,output] = fsolve(f,x0,options);

Finray1.plot_state(x_solve);

axis equal
axis([-80e-3 100e-3 0 160e-3]);


% 求出外力的合力
Fx=0;
Fy=0;
for i=1:N_tan
    Finray1.A_force_array(i).cal_pe;
    pka=Finray1.A_force_array(i).pe;

    PHI=pka(3)+Finray1.pA(3);
    Fx=Fx+tangent_F(i)*sin(PHI);
    Fy=Fy-tangent_F(i)*cos(PHI);
end

norm_F=norm([Fx;Fy]);
clc
disp(Fx)
disp(Fy)
disp(norm_F)


% 尝试在无梯度的情况下直接求解Finray与圆柱相切的问题

% 给定圆心位置x y 圆柱半径r

% 输入：相切位置L （以比例的形式给出，在0和1之间） 切点处的法向力F
% 输出：残差r1 r2
% 原理：使残差为0的L和F就是真实的情况

function r=myfunc(tangent_var,delta,center_x,center_y,a,b)

    N_tan=length(tangent_var)/3;
    
    tangent_t=zeros(N_tan,1);
    tangent_ratio=zeros(N_tan,1);
    tangent_F=zeros(N_tan,1);
    
    for i=1:N_tan
        tangent_t(i)=tangent_var(3*i-2);
        tangent_ratio(i)=tangent_var(3*i-1);
        tangent_F(i)=tangent_var(3*i);
    end

    
    xA=0-delta;
    yA=0;

    xB=73.24e-3-delta;
    yB=0;


    alpha_degree=90;
    beta_degree=116;
    psi_degree=beta_degree-alpha_degree;

    psi=psi_degree/180*pi;
    alpha=alpha_degree/180*pi;
    beta=beta_degree/180*pi;


    nA=50;
    nB=50;


    LA=(xB-xA)/tan(psi);
    LB=(xB-xA)/sin(psi);

    wid=28e-3;
    thi=0.2e-3;
    E=213e9;


    constraint_ratio=[];

    La_1=27e-3;
    La_2=54e-3;
    La_3=81e-3;
    La_4=108e-3;

    Lb_1=45.67e-3;
    Lb_2=78.69e-3;
    Lb_3=107.5e-3;
    Lb_4=132.27e-3;

    Lcon1=sqrt((LA-La_1)^2+(LB-Lb_1)^2-2*(LA-La_1)*(LB-Lb_1)*cos(psi));
    Lcon2=sqrt((LA-La_2)^2+(LB-Lb_2)^2-2*(LA-La_2)*(LB-Lb_2)*cos(psi));
    Lcon3=sqrt((LA-La_3)^2+(LB-Lb_3)^2-2*(LA-La_3)*(LB-Lb_3)*cos(psi));
    Lcon4=sqrt((LA-La_4)^2+(LB-Lb_4)^2-2*(LA-La_4)*(LB-Lb_4)*cos(psi));

%     constraint_ratio=[Lcon1,La_1/LA,Lb_1/LB;
%                       Lcon2,La_2/LA,Lb_2/LB;
%                       Lcon3,La_3/LA,Lb_3/LB;
%                       Lcon4,La_4/LA,Lb_4/LB];


    A_force_ratio=[tangent_F,tangent_ratio];


    B_force_ratio=[];


    % 根据上面设定的参数生成参数结构体
    finray_info=struct();

    finray_info.pA=[xA;yA;alpha];
    finray_info.pB=[xB;yB;beta];

    finray_info.LA=LA;
    finray_info.LB=LB;

    finray_info.nA=nA;
    finray_info.nB=nB;

    finray_info.psi=psi;

    finray_info.wid_A=wid;
    finray_info.wid_B=wid;

    finray_info.thi_A=thi;
    finray_info.thi_B=thi;

    finray_info.E_A=LA/135e-3*E;
    finray_info.E_B=LB/156e-3*E;

    finray_info.constraint_ratio=constraint_ratio;

    finray_info.A_force_ratio=A_force_ratio;

    finray_info.B_force_ratio=B_force_ratio;


    Finray1=finray_force(finray_info);




    x0=zeros(nA+nB+3+2*Finray1.constraint_number,1);

    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

    f=@(x) Finray1.cal_balance(x);
    [x_solve,fval,exitflag,output] = fsolve(f,x0,options);
    
    
    
    % 拿到切点的x y phi
    
    thetaA=x_solve(1:Finray1.nA);
    
    r=zeros(3*N_tan,1);
    
    for i=1:N_tan
        Finray1.A_force_array(i).theta=thetaA(1:Finray1.A_force_index(i,2));
        Finray1.A_force_array(i).cal_pe;

        pka=Finray1.A_force_array(i).pe;
        L_tail=Finray1.A_force_ratio(i,2)*Finray1.LA-(Finray1.A_force_array(i).Ltotal-Finray1.A_force_array(i).seg_length/2);
        pka(1)=pka(1)-Finray1.A_force_array(i).seg_length*cos(sum(Finray1.A_force_array(i).theta))/2+L_tail*cos(sum(Finray1.A_force_array(i).theta));
        pka(2)=pka(2)-Finray1.A_force_array(i).seg_length*sin(sum(Finray1.A_force_array(i).theta))/2+L_tail*sin(sum(Finray1.A_force_array(i).theta));

        PA(1)=pka(1)*cos(Finray1.pA(3))-pka(2)*sin(Finray1.pA(3))+Finray1.pA(1);
        PA(2)=pka(1)*sin(Finray1.pA(3))+pka(2)*cos(Finray1.pA(3))+Finray1.pA(2);
        PA(3)=pka(3)+Finray1.pA(3);
        
        r(3*i-2)=b*cos(tangent_t(i))*cos(PA(3))+a*sin(tangent_t(i))*sin(PA(3));
        r(3*i-1)=center_x+a*cos(tangent_t(i))-PA(1);
        r(3*i)=center_y+b*sin(tangent_t(i))-PA(2);
    end
    
    
end