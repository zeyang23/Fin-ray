% 验证demo17_4的计算结果

clear
clc

load('X_series_lsq.mat');

L0=1;

radius=0.1*L0;
center_x=0.2*L0;
center_y=0.5*L0;

delta=linspace(0.1*L0,0,21);


vid = VideoWriter('finray_contact_animation_lsq');
writerObj.FrameRate = 30;
open(vid);
    
for i=1:21
    rectangle('Position',[center_x-radius,center_y-radius,2*radius,2*radius],'Curvature',[1,1],'linewidth',1,'edgecolor','r')
    axis equal;
    
    hold on
    
    
    tangent_ratio=X_series_lsq(i,1);
    tangent_F=X_series_lsq(i,2);

    L0=1;

    xA=0+delta(i);
    yA=0;

    xB=0.35*L0+delta(i);
    yB=0*L0;

    psi_degree=20;
    alpha_degree=80;
    beta_degree=100;

    psi=psi_degree/180*pi;
    alpha=alpha_degree/180*pi;
    beta=beta_degree/180*pi;


    nA=50;
    nB=50;


    LA=1*L0;
    LB=1*L0;

    wid=5e-3;
    thi=1e-3;
    E=197*1e9;


    constraint_ratio=[];
    Lcon1=6/7*sqrt((xA-xB)^2+(yA-yB)^2);
    Lcon2=5/7*sqrt((xA-xB)^2+(yA-yB)^2);
    Lcon3=4/7*sqrt((xA-xB)^2+(yA-yB)^2);
    Lcon4=3/7*sqrt((xA-xB)^2+(yA-yB)^2);
    Lcon5=2/7*sqrt((xA-xB)^2+(yA-yB)^2);
    Lcon6=1/7*sqrt((xA-xB)^2+(yA-yB)^2);
    constraint_ratio=[Lcon1,1/7,1/7;
                  Lcon2,2/7,2/7;
                  Lcon3,3/7,3/7;
                  Lcon4,4/7,4/7;
                  Lcon5,5/7,5/7;
                  Lcon6,6/7,6/7];


    A_force_ratio=[];
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

    finray_info.E_A=E;
    finray_info.E_B=E;

    finray_info.constraint_ratio=constraint_ratio;

    finray_info.A_force_ratio=A_force_ratio;

    finray_info.B_force_ratio=B_force_ratio;


    Finray1=finray_force(finray_info);




    x0=zeros(nA+nB+3+2*Finray1.constraint_number,1);

    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

    g=@(x) Finray1.cal_balance(x);
    [x_solve,gval,exitflag_g,output2] = fsolve(g,x0,options);

    Finray1.plot_state(x_solve);
    
    axis([-0.1*L0,0.6*L0,0,L0])
    
    set(gcf,'doublebuffer','on');
    drawnow;    
    Frame = getframe;    
    writeVideo(vid,Frame);

    
    if i~=21
        clf;
    end
end

close(vid);