% 可变长度的柔性板，与圆柱相切
% 动画演示

clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;

delta=linspace(0.2*L0,0,41);

R1=planar_nR_flex(E,L0,wid,thi,n,[0;0;0]);

xA=0;
yA=0;

alpha_degree=80;
alpha=alpha_degree/180*pi;


radius=0.1*L0;
center_x=0.2*L0;
center_y=0.4*L0;




vid = VideoWriter('finray_contact_no_rigid');
writerObj.FrameRate = 30;
open(vid);

for i=1:length(delta)

    x0=zeros(n+2,1);
    x0(end-1)=L0;

    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

    f=@(x) myfunc(x,xA+delta(i),yA,alpha,radius,center_x,center_y,R1);
    [x_solve,~,~,~] = fsolve(f,x0,options);

    R1.cal_posall;
    plot_abs_pos(R1.pos_all,alpha,[xA+delta(i),yA]);
    hold on
    rectangle('Position',[center_x-radius,center_y-radius,2*radius,2*radius],'Curvature',[1,1],'linewidth',1,'edgecolor','r')
    axis equal;


    %%
    L0=1;

    R1.update;

    Temp=[cos(alpha) -sin(alpha);sin(alpha) cos(alpha)]*[R1.pe(1);R1.pe(2)]+[xA+delta(i);yA];

    xA2=Temp(1);
    yA2=Temp(2);

    alpha2=R1.pe(3)+alpha;
    alpha2_degree=alpha2/pi*180;

    xB=0.35*L0+delta(i);
    yB=0*L0;

    psi_degree=20;
    beta_degree=100;

    psi=psi_degree/180*pi;
    beta=beta_degree/180*pi;


    nA=50;
    nB=50;


    LA=1*L0-x_solve(end-1);
    LB=1*L0;

    wid=5e-3;
    thi=1e-3;
    E=197*1e9;


    constraint_ratio=[];
    A_force_ratio=[];
    B_force_ratio=[];



    % 根据上面设定的参数生成参数结构体
    finray_info=struct();

    finray_info.pA=[xA2;yA2;alpha2];
    finray_info.pB=[xB;yA;beta];

    finray_info.LA=LA;
    finray_info.LB=LB;

    finray_info.nA=nA;
    finray_info.nB=nB;

    finray_info.psi=psi;

    finray_info.wid_A=wid;
    finray_info.wid_B=wid;

    finray_info.thi_A=thi;
    finray_info.thi_B=thi;

    finray_info.E_A=E;
    finray_info.E_B=E;

    finray_info.constraint_ratio=constraint_ratio;

    finray_info.A_force_ratio=A_force_ratio;

    finray_info.B_force_ratio=B_force_ratio;


    Finray1=finray(finray_info);




    x0=zeros(nA+nB+3+2*Finray1.constraint_number,1);

    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

    f=@(x) Finray1.cal_balance(x);
    [x_solve,fval,exitflag,output] = fsolve(f,x0,options);

    hold on
    Finray1.plot_state(x_solve);
    
    
    
    axis([-0.1*L0,0.6*L0,0,L0])
    
    set(gcf,'doublebuffer','on');
    drawnow;    
    Frame = getframe;    
    writeVideo(vid,Frame);

    
    if i~=length(delta)
        clf;
    end

end

close(vid)

%%
function [r,J]=myfunc(x,xA,yA,alpha,radius,center_x,center_y,R1)
    
    thetaA=x(1:R1.n_seg);
    L_new=x(end-1);
    F_new=x(end);
    
    R1.theta=thetaA;
    R1.change_length(L_new);
   
    R1.F=[F_new*sin(R1.pe(3));-F_new*cos(R1.pe(3));0];
    
    R1.update;
    
    r=zeros(R1.n_seg+2,1);
    
    r(1:R1.n_seg)=R1.K_theta*R1.theta-transpose(R1.Jacobian)*R1.F;
    
    r(end-1:end)=[center_x-xA+radius*sin(R1.pe(3)+alpha)-(R1.pe(1)*cos(alpha)-R1.pe(2)*sin(alpha));
                  center_y-yA-radius*cos(R1.pe(3)+alpha)-(R1.pe(1)*sin(alpha)+R1.pe(2)*cos(alpha))];
              
    J=zeros(R1.n_seg+2);
    
    J(1:R1.n_seg,1:R1.n_seg)=R1.K_theta-R1.partial-transpose(R1.Jacobian)*[F_new*cos(R1.pe(3));F_new*sin(R1.pe(3));0]*ones(1,R1.n_seg);
    
    temp=transpose(R1.Jacobian);
    temp(:,3)=zeros(1,R1.n_seg);
    temp=1/L_new*temp;
    
    J(1:R1.n_seg,R1.n_seg+1)=-1/L_new*R1.K_theta*R1.theta-temp*R1.F;
    
    J(1:R1.n_seg,end)=-transpose(R1.Jacobian)*[sin(R1.pe(3));-cos(R1.pe(3));0];

    temp2=R1.Jacobian;
    J(end-1:end,1:R1.n_seg)=[radius*cos(R1.pe(3)+alpha)*ones(1,R1.n_seg);radius*sin(R1.pe(3)+alpha)*ones(1,R1.n_seg)]-[cos(alpha) -sin(alpha);sin(alpha) cos(alpha)]*temp2(1:2,:);
    
    J(end-1:end,end-1)=-[cos(alpha) -sin(alpha);sin(alpha) cos(alpha)]*1/L_new*[R1.pe(1);R1.pe(2)];
    
end