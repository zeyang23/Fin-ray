% 可变长度的柔性板，与圆柱相切

% 法向力为0

clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;


R1=planar_nR_flex(E,L0,wid,thi,n,[0;0;0]);

xA=0;
yA=0;

alpha_degree=80;
alpha=alpha_degree/180*pi;


radius=0.1*L0;
center_x=0.2*L0;
center_y=0.4*L0;


x0=zeros(n+2,1);
x0(end-1)=L0;

options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',true);

f=@(x) myfunc(x,xA,yA,alpha,radius,center_x,center_y,R1);
[x_solve,fval,exitflag,output] = fsolve(f,x0,options);

R1.cal_posall;
plot_abs_pos(R1.pos_all,alpha,[xA,yA]);
hold on
rectangle('Position',[center_x-radius,center_y-radius,2*radius,2*radius],'Curvature',[1,1],'linewidth',1,'edgecolor','r')
axis equal;




function [r,J]=myfunc(x,xA,yA,alpha,radius,center_x,center_y,R1)
    thetaA=x(1:R1.n_seg);
    L_new=x(end-1);
    F_new=x(end);
    
    R1.theta=thetaA;
    R1.change_length(L_new);
   
    R1.F=[F_new*cos(R1.pe(3));F_new*sin(R1.pe(3));0];
    
    R1.update;
    
    r=zeros(R1.n_seg+2,1);
    
    r(1:R1.n_seg)=R1.K_theta*R1.theta-transpose(R1.Jacobian)*R1.F;
    
    r(end-1:end)=[center_x-xA+radius*sin(R1.pe(3)+alpha)-(R1.pe(1)*cos(alpha)-R1.pe(2)*sin(alpha));
                  center_y-yA-radius*cos(R1.pe(3)+alpha)-(R1.pe(1)*sin(alpha)+R1.pe(2)*cos(alpha))];
              
    J=zeros(R1.n_seg+2);
    
    J(1:R1.n_seg,1:R1.n_seg)=R1.K_theta-R1.partial-transpose(R1.Jacobian)*[-F_new*sin(R1.pe(3));F_new*cos(R1.pe(3));0]*ones(1,R1.n_seg);
    
    temp=transpose(R1.Jacobian);
    temp(:,3)=zeros(1,R1.n_seg);
    temp=1/L_new*temp;
    
    J(1:R1.n_seg,R1.n_seg+1)=-1/L_new*R1.K_theta*R1.theta-temp*R1.F;
    
    J(1:R1.n_seg,end)=-transpose(R1.Jacobian)*[cos(R1.pe(3));sin(R1.pe(3));0];

    temp2=R1.Jacobian;
    J(end-1:end,1:R1.n_seg)=[radius*cos(R1.pe(3)+alpha)*ones(1,R1.n_seg);radius*sin(R1.pe(3)+alpha)*ones(1,R1.n_seg)]-[cos(alpha) -sin(alpha);sin(alpha) cos(alpha)]*temp2(1:2,:);
    
    J(end-1:end,end-1)=-[cos(alpha) -sin(alpha);sin(alpha) cos(alpha)]*1/L_new*[R1.pe(1);R1.pe(2)];
    
end