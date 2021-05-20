% 与圆的多点接触动画


%%
close all
clear
clc


radius=31.5e-3;
center_x=0-radius;
center_y=66.5e-3;

delta_total=36e-3;


% 相切点的个数
N_tan=3;


lb=zeros(2*N_tan,1);

ub=zeros(2*N_tan,1);
for i=1:N_tan
    ub(2*i-1)=1;
    ub(2*i)=10;
end

ratio_series=linspace(0.4,0.55,N_tan);

X0=zeros(2*N_tan,1);
for i=1:N_tan
    X0(2*i-1)=ratio_series(i);
    X0(2*i)=6;
end



N_total=37;
delta_series=linspace(0,delta_total,N_total);


X_series=[];
resnorm_series=[];
residual_series=[];

F_cal_series=[];


vid = VideoWriter('finray_circle');
writerObj.FrameRate = 30;
open(vid);

for i=1:N_total
    delta=delta_series(i);
    
    f=@(tangent_var) myfunc(tangent_var,delta,center_x,center_y,radius);

    [X,resnorm,residual,exitflag_lsq,output_lsq] = lsqnonlin(f,X0,lb,ub);
    
    
    X_series(:,i)=X;
    resnorm_series(:,i)=resnorm;
    residual_series(:,i)=residual;
    
    X0=X;
    X0(1)=X0(3)-0.2;
    X0(5)=X0(3)+0.2;


    %% 验证结果
    % 画出圆柱位置
    rectangle('Position',[center_x-radius,center_y-radius,2*radius,2*radius],'Curvature',[1,1],'linewidth',2,'edgecolor','r')
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

    constraint_ratio=[Lcon1,La_1/LA,Lb_1/LB;
                      Lcon2,La_2/LA,Lb_2/LB;
                      Lcon3,La_3/LA,Lb_3/LB;
                      Lcon4,La_4/LA,Lb_4/LB];

    tangent_ratio=zeros(N_tan,1);
    tangent_F=zeros(N_tan,1);

    for k=1:N_tan
        tangent_ratio(k)=X(2*k-1);
        tangent_F(k)=X(2*k);
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
    axis([-90e-3 100e-3 0 160e-3]);


    % 求出外力的合力
    Fx=0;
    Fy=0;
    for j=1:N_tan
        Finray1.A_force_array(j).cal_pe;
        pka=Finray1.A_force_array(j).pe;

        PHI=pka(3)+Finray1.pA(3);
        Fx=Fx+tangent_F(j)*sin(PHI);
        Fy=Fy-tangent_F(j)*cos(PHI);
    end

    norm_F=norm([Fx;Fy]);
%     clc
%     disp(Fx)
%     disp(Fy)
%     disp(norm_F)

    F_cal_series(:,i)=[Fx;Fy;norm_F];
    
    
    set(gcf,'doublebuffer','on');
    drawnow;    
    Frame = getframe;    
    writeVideo(vid,Frame);

    
    if i~=length(delta_series)
        clf;
    end
end

close(vid);

% 尝试在无梯度的情况下直接求解Finray与圆柱相切的问题

% 给定圆心位置x y 圆柱半径r

% 输入：相切位置L （以比例的形式给出，在0和1之间） 切点处的法向力F
% 输出：残差r1 r2
% 原理：使残差为0的L和F就是真实的情况

function r=myfunc(tangent_var,delta,center_x,center_y,radius)

    N_tan=length(tangent_var)/2;
    
    tangent_ratio=zeros(N_tan,1);
    tangent_F=zeros(N_tan,1);
    
    for i=1:N_tan
        tangent_ratio(i)=tangent_var(2*i-1);
        tangent_F(i)=tangent_var(2*i);
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

    constraint_ratio=[Lcon1,La_1/LA,Lb_1/LB;
                      Lcon2,La_2/LA,Lb_2/LB;
                      Lcon3,La_3/LA,Lb_3/LB;
                      Lcon4,La_4/LA,Lb_4/LB];


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
    
    r=zeros(2*N_tan,1);
    
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
        
        r(2*i-1)=center_x+radius*sin(PA(3))-PA(1);
        r(2*i)=center_y-radius*cos(PA(3))-PA(2);
    end
    
    
    
end