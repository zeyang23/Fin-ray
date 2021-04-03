% 可变长度的柔性板
% 与滑轮相切

% 21-04-03
% 第1部分（可变长度）和第2部分（不可变长度）都通过了Gradient Check

% 把梯度给fsolve比自己用牛顿法要强不少。。。
% 这波叫做把信息交给可靠的人 手动狗头
% 梯度在我手里我就只能玩牛顿法，梯度在fsolve手里可以被更好地利用


clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;


R1=planar_nR_flex(E,L0,wid,thi,n,[0;0;0]);

contact_degree=120;

pulley=[0.3*L0;0.3*L0;contact_degree/180*pi];


options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',true);


X0=zeros(n+3,1);
X0(end-2)=L0;

F=@(x) myfunc(x,pulley,R1);
[x_solve,fval,exitflag,output] = fsolve(F,X0,options);

R1.cal_posall;
R1.plot_all;


% 求解第2段

pos_end=[0.6*L0;0.6*L0;pi/2];

% 计算pos_end相对于pully的坐标

r_pulley=[cos(pulley(3)) sin(pulley(3));-sin(pulley(3)) cos(pulley(3))]*[-pulley(1)+pos_end(1);-pulley(2)+pos_end(2)];

phi_pulley=pos_end(3)-pulley(3);

pos_end_pulley=[r_pulley;phi_pulley];

n2=50;
R2=planar_nR(E,L0,wid,thi,n2,pos_end_pulley);

X0=zeros(n2+3,1);

F=@(x) myfunc2(x,R2);
[X_solve,Fval,Exitflag,Output] = fsolve(F,X0,options);

R2.cal_posall;
hold on
plot_abs_pos(R2.pos_all,pulley(3),[pulley(1),pulley(2)]);



% 画出滑轮
r=0.01*L0;
vec=r*[cos(pulley(3));sin(pulley(3));0];
vec1=rotz(90)*vec;
vec2=rotz(-90)*vec;
center1x=vec1(1)+pulley(1);
center1y=vec1(2)+pulley(2);
center2x=vec2(1)+pulley(1);
center2y=vec2(2)+pulley(2);

rectangle('Position',[center1x-r,center1y-r,2*r,2*r],'Curvature',[1,1],'linewidth',1,'edgecolor','r')
rectangle('Position',[center2x-r,center2y-r,2*r,2*r],'Curvature',[1,1],'linewidth',1,'edgecolor','r')



function [r,J]=myfunc2(x,R2)
    theta=x(1:R2.n_seg);
    F_new=x(end-2:end);
    
    R2.theta=theta;
    R2.F=F_new;
    
    R2.update;
    
    r=R2.cal_r;
    
    
    J=zeros(R2.n_seg+3);
    J(1:3,1:R2.n_seg)=R2.Jacobian;
    J(4:R2.n_seg+3,end-2:end)=-transpose(R2.Jacobian);

    temp=zeros(R2.n_seg);
    for i=1:R2.n_seg
        for j=1:R2.n_seg
            temp(i,j)=transpose(R2.F)*[R2.Hessian_xe(j,i);R2.Hessian_ye(j,i);0];
        end
    end

    J(4:end,1:R2.n_seg)=R2.K_theta-1*temp;
end


function [r,J]=myfunc(x,pulley,R1)
    



    thetaA=x(1:R1.n_seg);
    L_new=x(end-2);
    F_new=x(end-1);
    m_new=x(end);
    
    
    R1.theta=thetaA;
    R1.change_length(L_new);
   
    R1.F=[F_new*sin(R1.pe(3));-F_new*cos(R1.pe(3));m_new];
    
    R1.update;
    
    
    r=zeros(R1.n_seg+3,1);
    
    r(1:R1.n_seg)=R1.K_theta*R1.theta-transpose(R1.Jacobian)*R1.F;
    
    r(end-2:end)=R1.pe-pulley;
    
    
    J=zeros(R1.n_seg+3);
    
    J(1:R1.n_seg,1:R1.n_seg)=R1.K_theta-R1.partial-transpose(R1.Jacobian)*[F_new*cos(R1.pe(3));F_new*sin(R1.pe(3));0]*ones(1,R1.n_seg);
    
    temp=transpose(R1.Jacobian);
    temp(:,3)=zeros(1,R1.n_seg);
    temp=1/L_new*temp;
    
    J(1:R1.n_seg,R1.n_seg+1)=-1/L_new*R1.K_theta*R1.theta-temp*R1.F;
    
    J(1:R1.n_seg,end-1)=-transpose(R1.Jacobian)*[sin(R1.pe(3));-cos(R1.pe(3));0];
    
    J(1:R1.n_seg,end)=-transpose(R1.Jacobian)*[0;0;1];
    

    J(end-2:end,1:R1.n_seg)=R1.Jacobian;
    
    J(end-2:end,end-2)=[1/L_new*R1.pe(1);1/L_new*R1.pe(2);0];
    
end