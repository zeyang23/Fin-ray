% 无刚性约束的finray
% 静摩擦
% 临界条件：法向力为0 待验证


%% 初始状态
clear
clc

wid=14e-3;
thi=0.3e-3;
E=197*1e9;

Iz=1/12*thi.^3*wid;


L0=0.25;

LA=L0;
LB=L0;


IA=Iz;
IB=Iz;


half_psi_degree=10;
psi_degree=20;


alpha_degree=90-half_psi_degree;
beta_degree=180-alpha_degree;

psi=psi_degree/180*pi;
alpha=alpha_degree/180*pi;
beta=beta_degree/180*pi;


xA=0.5;
yA=0;

xB=0.5+2*L0*sin(psi/2);
yB=0;

nA=50;
nB=50;

problem_info_init.wid=wid;
problem_info_init.thi=thi;
problem_info_init.nA=nA;
problem_info_init.nB=nB;
problem_info_init.LA=LA;
problem_info_init.LB=LB;
problem_info_init.E=E;
problem_info_init.alpha=alpha;
problem_info_init.beta=beta;
problem_info_init.psi=psi;
problem_info_init.xA=xA;
problem_info_init.yA=yA;
problem_info_init.xB=xB;
problem_info_init.yB=yB;

solve_finray(problem_info_init,1);


X0=0.5+L0/2*sin(psi/2);
Y0=L0/2*cos(psi/2);
Theta_0=-psi/2;

radius=0.05;
Xc=X0-radius*cos(Theta_0);
Yc=Y0-radius*sin(Theta_0);

L0_up=L0/2;
L0_down=L0/2;

hold on
rectangle('Position',[Xc-radius,Yc-radius,2*radius,2*radius],'Curvature',[1,1],'linewidth',1,'edgecolor','r')

%%
DELTA=0.03;
N=40;
DELTA_series=linspace(0,DELTA,N+1);

Theta0_up=Theta_0;
Theta0_down=Theta_0;

theta_up_series=[];
res_up_series=[];

theta_down_series=[];
res_down_series=[];

for i=1:length(DELTA_series)
    delta=DELTA_series(i);
    f1=@(theta_up) solve_up(theta_up,delta,Xc,Yc,radius,Theta0_up,L0_up,problem_info_init,0);

    lb_1=Theta0_up;
    ub_1=pi;
    [theta_up_solve,resnorm,residual,exitflag_lsq,output_lsq] = lsqnonlin(f1,Theta0_up,lb_1,ub_1);
%     [theta_up_solve,fval_fsolve,exitflag_fsolve,output_fsolve] = fsolve(f1,Theta_0);

    theta_up_series(i)=theta_up_solve;
    res_up_series(i)=resnorm;

    X_up=Xc+radius*cos(theta_up_solve);
    Y_up=Yc+radius*sin(theta_up_solve);

    problem_info=problem_info_init;
    problem_info.xA=X_up;
    problem_info.yA=Y_up;
    problem_info.alpha=theta_up_solve+pi/2;

    problem_info.xB=problem_info.xB-delta;

    problem_info.LA=L0_up-radius*(theta_up_solve-Theta0_up);
    
    

    x=solve_finray(problem_info,1);
    clear problem_info
    
    L0_up=L0_up-radius*(theta_up_solve-Theta0_up);
    Theta0_up=theta_up_solve;
    
    
    
%     f2=@(theta_down) solve_down(theta_down,DELTA,Xc,Yc,radius,Theta0_down,L0_down,problem_info_init,0);
%     
%     lb_2=-pi;
%     ub_2=Theta0_down;
%     [theta_down_solve,resnorm,residual,exitflag_lsq,output_lsq] = lsqnonlin(f2,Theta0_down,lb_2,ub_2);
% %     [theta_up_solve,fval_fsolve,exitflag_fsolve,output_fsolve] = fsolve(f2,Theta0_down);
% %     resnorm=norm(fval_fsolve);
%     
%     theta_down_series(i)=theta_down_solve;
%     res_down_series(i)=resnorm;
% 
%     X_down=Xc+radius*cos(theta_down_solve);
%     Y_down=Yc+radius*sin(theta_down_solve);
%     
%     problem_info=problem_info_init;
%     problem_info.x_des=X_down;
%     problem_info.y_des=Y_down;
%     problem_info.phi_des=theta_down_solve+pi/2;
%     
%     problem_info.xA=problem_info.xA-delta;
%     
%     problem_info.LA=L0_down-radius*(Theta0_down-theta_down_solve);
%     
%     x=solve_rod(problem_info,1);
%     clear problem_info
%     
%     L0_down=L0_down-radius*(Theta0_down-theta_down_solve);
%     Theta0_down=theta_down_solve;
    
    
    

    
    rectangle('Position',[Xc-radius,Yc-radius,2*radius,2*radius],'Curvature',[1,1],'linewidth',1,'edgecolor','r')
    pause(0.02);
    
    if i~=length(DELTA_series)
        clf;
    end
    
end


%%
function r=solve_down(theta_down,DELTA,Xc,Yc,radius,Theta0_down,L0_down,problem_info_init,plot_flag)
    X_down=Xc+radius*cos(theta_down);
    Y_down=Yc+radius*sin(theta_down);
    
    problem_info=problem_info_init;
    problem_info.x_des=X_down;
    problem_info.y_des=Y_down;
    problem_info.phi_des=theta_down+pi/2;
    
    problem_info.xA=problem_info.xA-DELTA;
    
    problem_info.LA=L0_down-radius*(Theta0_down-theta_down);
    
    x=solve_rod(problem_info,plot_flag);
    
    alpha=problem_info.alpha;
    
    r=[cos(theta_down) sin(theta_down)]*[cos(alpha) -sin(alpha);sin(alpha) cos(alpha)]*[x(end-2);x(end-1)];
end

function r=solve_up(theta_up,DELTA,Xc,Yc,radius,Theta_0,L0_up,problem_info_init,plot_flag)
    X_up=Xc+radius*cos(theta_up);
    Y_up=Yc+radius*sin(theta_up);
    
    problem_info=problem_info_init;
    problem_info.xA=X_up;
    problem_info.yA=Y_up;
    problem_info.alpha=theta_up+pi/2;
    
    problem_info.xB=problem_info.xB-DELTA;
    
    problem_info.LA=L0_up-radius*(theta_up-Theta_0);
    
    x=solve_finray(problem_info,plot_flag);
    
    beta=problem_info.beta;
    
    r=[cos(theta_up) sin(theta_up)]*[cos(beta) -sin(beta);sin(beta) cos(beta)]*[x(end-2);x(end-1)];
    
end


function x=solve_rod(problem_info,plot_flag)
    L0=problem_info.LA;
    
    E=problem_info.E;
    
    wid=problem_info.wid;
    thi=problem_info.thi;
    nA=problem_info.nA;
    
    xA=problem_info.xA;
    yA=problem_info.yA;
    alpha=problem_info.alpha;
    
    pA=[xA;yA;alpha];

    x_des=problem_info.x_des;
    y_des=problem_info.y_des;
    phi_des=problem_info.phi_des;
    
    pdes=[x_des;y_des;phi_des];
    
    R1=planar_nR(E,L0,wid,thi,nA,pdes);
    
    f=@(x) cal_balance(x,R1,pA,pdes);
    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');
    x0=zeros(R1.n_seg+3,1);
    [x_solve,fval,exitflag,output] = fsolve(f,x0,options);
    
    x=x_solve;
    
    if exitflag<=0
        error('fail')
    end

    R1.theta=x_solve(1:R1.n_seg);
    R1.F=x_solve(R1.n_seg+1:R1.n_seg+3);

    R1.update;
    R1.cal_posall;
    
    if plot_flag ==1

        plot_abs_pos(R1.pos_all,alpha,[xA,yA]);
    end
end

function x=solve_finray(problem_info,plot_flag)
    wid=problem_info.wid;
    thi=problem_info.thi;
    nA=problem_info.nA;
    nB=problem_info.nB;
    
    LA=problem_info.LA;
    LB=problem_info.LB;
    
    E=problem_info.E;
    
    alpha=problem_info.alpha;
    beta=problem_info.beta;
    psi=problem_info.psi;
    
    xA=problem_info.xA;
    yA=problem_info.yA;
    xB=problem_info.xB;
    yB=problem_info.yB;
    
    
    alpha_degree=alpha/pi*180;
    beta_degree=beta/pi*180;


    A=rotz(-alpha_degree)*rotz(beta_degree);
    b=rotz(-alpha_degree)*[xB-xA;yB-yA;beta-alpha-psi];
    B=-rotz(beta_degree-alpha_degree);


    pdes=[0;0;0];
    RA=planar_nR(E,LA,wid,thi,nA,pdes);
    RB=planar_nR(E,LB,wid,thi,nB,pdes);


    % 牛顿法求解2N+3方程组
    x0=zeros(nA+nB+3,1);

    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

    f=@(x) myfun(x,nA,nB,A,B,b,RA,RB);
    [x,fval,exitflag,output] = fsolve(f,x0,options);
    
    if exitflag<=0
        error('fail')
    end

    RA.cal_posall;
    RB.cal_posall;

    A_pos_all=RA.pos_all;
    B_pos_all=RB.pos_all;
    
    if plot_flag==1
        
        A_abs_pos_all=plot_abs_pos(A_pos_all,alpha,[xA,yA]);
        hold on
        B_abs_pos_all=plot_abs_pos(B_pos_all,beta,[xB,yB]);
    end
end


function [r,J]=myfun(x,nA,nB,A,B,b,RA,RB)
    
    thetaA=zeros(nA,1);
    thetaB=zeros(nB,1);
    FB=zeros(3,1);
    
    thetaA=x(1:nA);
    thetaB=x(nA+1:nA+nB);
    FB=x(end-2:end);
    
    FA=B*FB;
    
    RA.theta=thetaA;
    RA.F=FA;
    RB.theta=thetaB;
    RB.F=FB;
    
    RA.update;
    RB.update;
    
    r1=RA.pe-A*RB.pe-b;
    r2=RA.K_theta*thetaA-transpose(RA.Jacobian)*FA;
    r3=RB.K_theta*thetaB-transpose(RB.Jacobian)*FB;
    r=[r1;r2;r3];
    
    if nargout > 1
        J=zeros(nA+nB+3);

        J(1:3,1:nA)=RA.Jacobian;
        J(1:3,nA+1:nA+nB)=-A*RB.Jacobian;
        J(4:nA+3,1:nA)=RA.K_theta-1*RA.partial;
        J(4:nA+3,end-2:end)=-transpose(RA.Jacobian)*B;
        J(nA+4:end,nA+1:nA+nB)=RB.K_theta-1*RB.partial;
        J(nA+4:end,end-2:end)=-transpose(RB.Jacobian);
    end
end


function [r,J]=cal_balance(x,Rod,pA,pdes)
    alpha=pA(3);
    
    alpha_degree=alpha/pi*180;
    
    Rod.theta=x(1:Rod.n_seg);
    Rod.F=x(Rod.n_seg+1:Rod.n_seg+3);
    
    Rod.update;
    
    r=zeros(size(x));
    J=zeros(length(x));
    
    r(1:3)=rotz(alpha_degree)*Rod.pe+pA-pdes;
    r(4:Rod.n_seg+3)=Rod.K_theta*Rod.theta-transpose(Rod.Jacobian)*Rod.F;
    
    J=Rod.cal_rdot;
    J(1:3,1:Rod.n_seg)=rotz(alpha_degree)*Rod.Jacobian;
    J(4:end,1:Rod.n_seg)=Rod.K_theta-Rod.partial;
    J(4:end,end-2:end)=-transpose(Rod.Jacobian);
    
end