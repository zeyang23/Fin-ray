% 根据保存的结果来提取刚性支撑连接点的坐标，用以论文画图

close all;

clear all;
clc

load('X_series_r22.5_delta30.mat')


delta_total=30e-3;
N_total=size(X_series,2);
N_tan=3;

delta_series=linspace(0,delta_total,N_total);



i=31;
    
delta=delta_series(i);
X=X_series(:,i);


X0=X;
X0(1)=X0(3)-0.2;
X0(5)=X0(3)+0.2;




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
axis([-80e-3 100e-3 0 160e-3]);

xlabel('x/m')
ylabel('y/m')


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

rigid_abs_pos=Finray1.get_rigid_pos(x_solve);
rigid_abs_pos(:,1)=rigid_abs_pos(:,1)-xA;
rigid_abs_pos=1000*rigid_abs_pos;


print('r20_delta30','-djpeg','-r600')
save('r20_delta30_cal.mat','rigid_abs_pos')