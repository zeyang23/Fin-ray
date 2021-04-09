% 逐渐增大的外力
% 相应的变形动画

clear
clc

L0=1;

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

Finray1=finray(finray_info);

N=100;
F=5;
Fseries=linspace(0,F,N);
Fseries=[Fseries,flip(Fseries)];
x0=zeros(nA+nB+3+2*Finray1.constraint_number,1);


vid = VideoWriter('finray_animation');
writerObj.FrameRate = 60;
open(vid);


for i=1:2*N
    Fi=Fseries(i);
    A_force_ratio=[Fi,1/2];
    finray_info.A_force_ratio=A_force_ratio;
    
    Finray1=finray(finray_info);

    

    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

    f=@(x) Finray1.cal_balance(x);
    [x_solve,fval,exitflag,output] = fsolve(f,x0,options);

    Finray1.plot_state(x_solve);
    axis([-0.2*L0,0.6*L0,0,L0]);
    
    set(gcf,'doublebuffer','on');
%     pause(0.01);
    drawnow;    
    Frame = getframe;    
    writeVideo(vid,Frame);
    

    clf;

    
    x0=x_solve;
end


for i=1:2*N
    Fi=Fseries(i);
    B_force_ratio=[Fi,1/2];
    finray_info.B_force_ratio=B_force_ratio;
    
    Finray1=finray(finray_info);

    

    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

    f=@(x) Finray1.cal_balance(x);
    [x_solve,fval,exitflag,output] = fsolve(f,x0,options);

    Finray1.plot_state(x_solve);
    axis([-0.2*L0,0.6*L0,0,L0]);
    
    set(gcf,'doublebuffer','on');
%     pause(0.01);
    drawnow;    
    Frame = getframe;    
    writeVideo(vid,Frame);
    
    if i ~= 2*N
        clf;
    end
    
    x0=x_solve;
end

close(vid);