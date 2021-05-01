% 5月1日实验
% Finray

%% point load animation
clear
clc

xA=0;
yA=0;

xB=73.24e-3;
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
E=197e9;


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

finray_info.A_force_ratio=[];
finray_info.B_force_ratio=[];

Finray1=finray(finray_info);

N=100;
FA=8;
Fseries_A=linspace(0,FA,N);
Fseries_A=[Fseries_A,flip(Fseries_A)];

FB=7;
Fseries_B=linspace(0,FB,N);
Fseries_B=[Fseries_B,flip(Fseries_B)];
x0=zeros(nA+nB+3+2*Finray1.constraint_number,1);


vid = VideoWriter('test_05_01_finray_point_load_animation');
writerObj.FrameRate = 60;
open(vid);


for i=1:2*N
    Fi=Fseries_A(i);
    A_force_ratio=[Fi,(0.5*La_2+0.5*La_3)/LA];
    finray_info.A_force_ratio=A_force_ratio;
    
    Finray1=finray(finray_info);

    

    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

    f=@(x) Finray1.cal_balance(x);
    [x_solve,fval,exitflag,output] = fsolve(f,x0,options);

    Finray1.plot_state(x_solve);
    axis equal
    
    axis([-50e-3 120e-3 0 160e-3]);
    
    set(gcf,'doublebuffer','on');
%     pause(0.01);
    drawnow;    
    Frame = getframe;    
    writeVideo(vid,Frame);
    

    clf;

    
    x0=x_solve;
end


for i=1:2*N
    Fi=Fseries_B(i);
    B_force_ratio=[Fi,(0.5*Lb_2+0.5*Lb_3)/LB];
    finray_info.B_force_ratio=B_force_ratio;
    
    Finray1=finray(finray_info);

    

    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off');

    f=@(x) Finray1.cal_balance(x);
    [x_solve,fval,exitflag,output] = fsolve(f,x0,options);

    Finray1.plot_state(x_solve);
    axis equal
    
    axis([-50e-3 120e-3 0 160e-3]);
    
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