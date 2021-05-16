% 使用指数坐标法求解有一根刚性支撑的Fin-ray
% 手写牛顿法

% 21-05-16
% 效果不错


clear
clc

L0=1;

LA=0.8*L0;
nA=50;

wid_A=5e-3;
thi_A=1e-3;
E_A=197*1e9;

Iz_A=1/12*thi_A.^3*wid_A;
seg_length_A=LA/nA;
K_A=diag((E_A*Iz_A/seg_length_A)*ones(nA,1));


LB=L0;
nB=50;

wid_B=5e-3;
thi_B=1e-3;
E_B=197*1e9;

Iz_B=1/12*thi_B.^3*wid_B;
seg_length_B=LB/nB;
K_B=diag((E_B*Iz_B/seg_length_B)*ones(nB,1));



alpha_degree=80;
alpha=alpha_degree/180*pi;

beta_degree=100;
beta=beta_degree/180*pi;

psi=beta-alpha;


xA=0;
yA=0;
xB=0.35*L0;
yB=0;


L_left=0.5*LA;
L_right=0.5*LB;

Lcon=0.5*sqrt((xA-xB)^2+(yA-yB)^2);

N1=fix((L_left-seg_length_A/2)/seg_length_A)+1;
N2=fix((L_right-seg_length_B/2)/seg_length_B)+1;


Slist_A=zeros(6,nA);
Qlist_A=zeros(3,nA);

for i=1:nA
    q=[xA;yA;0]+seg_length_A/2*[cos(alpha);sin(alpha);0]+(i-1)*seg_length_A*[cos(alpha);sin(alpha);0];
    Qlist_A(:,i)=q;
    Slist_A(:,i)=[0;0;1;-cross([0;0;1],q)];
end


Slist_B=zeros(6,nB);
Qlist_B=zeros(3,nB);

for i=1:nB
    q=[xB;yB;0]+seg_length_B/2*[cos(beta);sin(beta);0]+(i-1)*seg_length_B*[cos(beta);sin(beta);0];
    Qlist_B(:,i)=q;
    Slist_B(:,i)=[0;0;1;-cross([0;0;1],q)];
end

g0_A=[rotz(alpha_degree),[xA;yA;0]+LA*[cos(alpha);sin(alpha);0];[0,0,0,1]];
g0_B=[rotz(beta_degree),[xB;yB;0]+LB*[cos(beta);sin(beta);0];[0,0,0,1]];

g0_F1=[rotz(alpha_degree),[xA;yA;0]+L_left*[cos(alpha);sin(alpha);0];[0,0,0,1]];
g0_F2=[rotz(beta_degree),[xB;yB;0]+L_right*[cos(beta);sin(beta);0];[0,0,0,1]];


problem_info=struct;
problem_info.LA=LA;
problem_info.nA=nA;
problem_info.seg_length_A=seg_length_A;
problem_info.K_A=K_A;

problem_info.LB=LB;
problem_info.nB=nB;
problem_info.seg_length_B=seg_length_B;
problem_info.K_B=K_B;

problem_info.alpha=alpha;
problem_info.beta=beta;
problem_info.psi=psi;

problem_info.Slist_A=Slist_A;
problem_info.Slist_B=Slist_B;

problem_info.g0_A=g0_A;
problem_info.g0_B=g0_B;
problem_info.g0_F1=g0_F1;
problem_info.g0_F2=g0_F2;

problem_info.Lcon=Lcon;
problem_info.N1=N1;
problem_info.N2=N2;


x=zeros(nA+nB+8,1);
TOL=1e-8;

k=1;
while(1)
    if k>500
        error("can not converge")
    end
    [b,J]=cal_balance(x,problem_info);
    
    delta=-pinv(J)*b;

    % 更新conv坐标系下的参数
    x=x+delta;

    k=k+1;
%         if(norm(x(:,k)-x(:,k-1))<TOL)
    if(norm(b)<TOL)
        fprintf('Newton Method converge: iteration = %d\n',k-1)
        fprintf('norm(e) = %E\n',norm(b))
        break;
    end
end



% 画出求解结果
theta_A=x(1:nA);
theta_B=x(nA+1:nA+nB);

g_A=FKinSpace(g0_A,Slist_A,theta_A);
g_B=FKinSpace(g0_B,Slist_B,theta_B);

g_F1=FKinSpace(g0_F1,Slist_A(:,1:N1),theta_A(1:N1));
g_F2=FKinSpace(g0_F2,Slist_B(:,1:N2),theta_B(1:N2));

r1=[eye(3),zeros(3,1)]*g_F1*[0;0;0;1];
r2=[eye(3),zeros(3,1)]*g_F2*[0;0;0;1];


plot(xA,yA,'o')
hold on
plot(xB,yB,'o')

plot(g_A(1,4),g_A(2,4),'o')
plot(g_B(1,4),g_B(2,4),'o')

pos_all_A=cal_axes_pos(theta_A,Qlist_A,Slist_A);
pos_all_B=cal_axes_pos(theta_B,Qlist_B,Slist_B);

plot(pos_all_A(1,:),pos_all_A(2,:))
plot(pos_all_B(1,:),pos_all_B(2,:))

plot([r1(1) r2(1)],[r1(2) r2(2)])

axis equal

function [res,J]=cal_balance(x,problem_info)
    nA=problem_info.nA;
    nB=problem_info.nB;
    
    N1=problem_info.N1;
    N2=problem_info.N2;
    
    Lcon=problem_info.Lcon;
    
    
    theta_A=x(1:nA);
    theta_B=x(nA+1:nA+nB);
    Fe=x(nA+nB+1:nA+nB+6);
    fcon=x(nA+nB+7);
    gamma=x(nA+nB+8);
    
    
    g_A=FKinSpace(problem_info.g0_A,problem_info.Slist_A,theta_A);
    g_B=FKinSpace(problem_info.g0_B,problem_info.Slist_B,theta_B);
    
    Rz_psi=[rotz(problem_info.psi/pi*180),[0;0;0];[0 0 0 1]];
    
    res1=se3ToVec(logm(g_A*Rz_psi*TransInv(g_B)));
    
    Js_A=JacobianSpace(problem_info.Slist_A,theta_A);
    Js_B=JacobianSpace(problem_info.Slist_B,theta_B);
    
    Jb_B=Adjoint(TransInv(g_B))*Js_B;
    
    J_F1=zeros(size(Js_A));
    J_F1(:,1:N1)=Js_A(:,1:N1);
    J_F2=zeros(size(Js_B));
    J_F2(:,1:N2)=Js_B(:,1:N2);
    
    f1=[-fcon*cos(gamma);-fcon*sin(gamma);0];
    f2=-f1;
    
    g_F1=FKinSpace(problem_info.g0_F1,problem_info.Slist_A(:,1:N1),theta_A(1:N1));
    g_F2=FKinSpace(problem_info.g0_F2,problem_info.Slist_B(:,1:N2),theta_B(1:N2));
    
    r1=[eye(3),zeros(3,1)]*g_F1*[0;0;0;1];
    r2=[eye(3),zeros(3,1)]*g_F2*[0;0;0;1];
    
    F1=[cross(r1,f1);f1];
    F2=[cross(r2,f2);f2];
    
    res2=problem_info.K_A*theta_A-transpose(Js_A)*Fe-transpose(J_F1)*F1;
    res3=problem_info.K_B*theta_B+transpose(Js_B)*Fe-transpose(J_F2)*F2;
    res4=[1 0 0;0 1 0]*(r2-r1)-Lcon*[cos(gamma);sin(gamma)];
    
    res=[res1;res2;res3;res4];
    
    J=zeros(length(x));
    
    J(1:6,1:nA)=Js_A;
    J(1:6,nA+1:nA+nB)=-Adjoint(g_A*Rz_psi)*Jb_B;
    
    
    
    partial_r1_thetaA=zeros(3,nA);
    for i=1:N1
        partial_r1_thetaA(:,i)=[eye(3),zeros(3,1)]*VecTose3(Js_A(:,i))*g_F1*[0;0;0;1];
    end
    
    partial_r2_thetaB=zeros(3,nB);
    for i=1:N2
        partial_r2_thetaB(:,i)=[eye(3),zeros(3,1)]*VecTose3(Js_B(:,i))*g_F2*[0;0;0;1];
    end
    
    KF_a=zeros(nA);
    for i=1:N1
        KF_a(:,i)=transpose(J_F1)*[cross(partial_r1_thetaA(:,i),f1);[0;0;0]];
    end
    
    KF_b=zeros(nB);
    for i=1:N2
        KF_b(:,i)=transpose(J_F2)*[cross(partial_r2_thetaB(:,i),f2);[0;0;0]];
    end

    
    J(7:nA+6,1:nA)=problem_info.K_A-get_partial(Js_A,Fe)-get_partial(J_F1,F1)-KF_a;
    J(7:nA+6,nA+nB+1:nA+nB+6)=-transpose(Js_A);
    J(7:nA+6,nA+nB+7)=-transpose(J_F1)*[cross(r1,[-cos(gamma);-sin(gamma);0]);[-cos(gamma);-sin(gamma);0]];
    J(7:nA+6,nA+nB+8)=-transpose(J_F1)*[cross(r1,[fcon*sin(gamma);-fcon*cos(gamma);0]);[fcon*sin(gamma);-fcon*cos(gamma);0]];
    
    
    J(nA+7:nA+nB+6,nA+1:nA+nB)=problem_info.K_B+get_partial(Js_B,Fe)-get_partial(J_F2,F2)-KF_b;
    J(nA+7:nA+nB+6,nA+nB+1:nA+nB+6)=transpose(Js_B);
    J(nA+7:nA+nB+6,nA+nB+7)=-transpose(J_F2)*[cross(r2,[cos(gamma);sin(gamma);0]);[cos(gamma);sin(gamma);0]];
    J(nA+7:nA+nB+6,nA+nB+8)=-transpose(J_F2)*[cross(r2,[-fcon*sin(gamma);fcon*cos(gamma);0]);[-fcon*sin(gamma);fcon*cos(gamma);0]];
    
    
    J(nA+nB+7:nA+nB+8,1:nA)=[1 0 0;0 1 0]*(-partial_r1_thetaA);
    J(nA+nB+7:nA+nB+8,nA+1:nA+nB)=[1 0 0;0 1 0]*partial_r2_thetaB;
    J(nA+nB+7:nA+nB+8,nA+nB+8)=Lcon*[sin(gamma);-cos(gamma)];
end

function partial=get_partial(J,F)
    n=size(J,2);
    partial=zeros(n);
    
    for i=2:n
        for j=1:i-1
            partial(i,j)=transpose(F)*ad(J(:,j))*J(:,i);
        end
    end
end