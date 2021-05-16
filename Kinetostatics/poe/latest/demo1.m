% 指数坐标法求解无刚性支撑的Finray
% 使用fsolve

% 21-05-16
% 使用fsolve求解效果很差劲

clear
clc

L0=1;

LA=L0;
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



alpha_degree=60;
alpha=alpha_degree/180*pi;

beta_degree=120;
beta=beta_degree/180*pi;

psi=beta-alpha;


xA=0;
yA=0;
xB=L0;
yB=0;



Slist_A=zeros(6,nA);

for i=1:nA
    q=[xA;yA;0]+seg_length_A/2*[cos(alpha);sin(alpha);0]+(i-1)*seg_length_A*[cos(alpha);sin(alpha);0];
    Slist_A(:,i)=[0;0;1;-cross([0;0;1],q)];
end


Slist_B=zeros(6,nB);

for i=1:nB
    q=[xB;yB;0]+seg_length_B/2*[cos(beta);sin(beta);0]+(i-1)*seg_length_B*[cos(beta);sin(beta);0];
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


x0=zeros(nA+nB+6,1);

f=@(x) cal_balance(x,problem_info);

options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',true);
[x_solve,fval,exitflag,output] = fsolve(f,x0,options);

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