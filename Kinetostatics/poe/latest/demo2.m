% 指数坐标法求解无刚性支撑的Finray
% 使用手写的牛顿法

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


x=zeros(nA+nB+6,1);
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



plot(xA,yA,'o')
hold on
plot(xB,yB,'o')

plot(g_A(1,4),g_A(2,4),'o')
plot(g_B(1,4),g_B(2,4),'o')

pos_all_A=cal_axes_pos(theta_A,Qlist_A,Slist_A);
pos_all_B=cal_axes_pos(theta_B,Qlist_B,Slist_B);

plot(pos_all_A(1,:),pos_all_A(2,:))
plot(pos_all_B(1,:),pos_all_B(2,:))

axis equal

function [r,J]=cal_balance(x,problem_info)
    nA=problem_info.nA;
    nB=problem_info.nB;
    
    theta_A=x(1:nA);
    theta_B=x(nA+1:nA+nB);
    Fe=x(nA+nB+1:nA+nB+6);
    
    g_A=FKinSpace(problem_info.g0_A,problem_info.Slist_A,theta_A);
    g_B=FKinSpace(problem_info.g0_B,problem_info.Slist_B,theta_B);
    
    Rz_psi=[rotz(problem_info.psi/pi*180),[0;0;0];[0 0 0 1]];
    
    r1=se3ToVec(logm(g_A*Rz_psi*TransInv(g_B)));
    
    Js_A=JacobianSpace(problem_info.Slist_A,theta_A);
    Js_B=JacobianSpace(problem_info.Slist_B,theta_B);
    
    Jb_B=Adjoint(TransInv(g_B))*Js_B;
    
    r2=problem_info.K_A*theta_A-transpose(Js_A)*Fe;
    r3=problem_info.K_B*theta_B+transpose(Js_B)*Fe;
    
    r=[r1;r2;r3];
    
    J=zeros(length(x));
    
    J(1:6,1:nA)=Js_A;
    J(1:6,nA+1:nA+nB)=-Adjoint(g_A*Rz_psi)*Jb_B;

    
    J(7:nA+6,1:nA)=problem_info.K_A-get_partial(Js_A,Fe);
    J(7:nA+6,nA+nB+1:nA+nB+6)=-transpose(Js_A);
    
    J(nA+7:nA+nB+6,nA+1:nA+nB)=problem_info.K_B+get_partial(Js_B,Fe);
    J(nA+7:nA+nB+6,nA+nB+1:nA+nB+6)=transpose(Js_B);
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