%% 2-a)
clear all;clc;addpath(genpath('.'));     %标准开头
a=[0 0.5 0 -pi/2 0];b=[0 0 0.5 0 0];
c=[0 0 0 -pi/2 0];d=[0 0.6 0 pi/2 0];            %设置6杆参数,长度单位化
e=[0 0 0 -pi/2 0];f=[0 0.5 0 0 0];
L1=Link(a);L2=Link(b);L3=Link(c);
L4=Link(d);L5=Link(e);L6=Link(f);
rob=SerialLink([L1,L2,L3,L4,L5,L6],'offset',[0 0 -pi/2 0 0 0]); 
g0=[[0 0 1;0 -1 0;1 0 0],[1.6;0;0.5];0,0,0,1];                    %初始姿态矩阵
w_all=[[0;0;1],[0;1;0],[0;1;0],[1;0;0],[0;1;0],[1;0;0]];    %所有的w向量构成矩阵
q_all=[[0;0;0],[0;0;0.5],[0.5;0;0.5],[0;0;0.5],[1.1;0;0.5],[0;0;0.5]];
theta_all=pi/20*[1;2;3;4;5;6];   %所有的q向量构成矩阵
[g_st,xi_all]=exp_fkine(w_all,q_all,g0,theta_all);   %指数正向运算，输入w_all，q_all，g0，theta_all
g_rtb=rob.fkine(theta_all);

rmpath(genpath('.'));                     %标准结局

 %% functions
function [g_st,xi_all] = exp_fkine(w_all,q_all,g0,theta_all)%正向运动学函数
xi_all=hight(w_all,q_all);                          %将R^6→se(3)
gi_mult=mult(xi_all,theta_all);                     %将se(3)→SE(3)并累乘
g_st=gi_mult*g0;                                    %初值矩阵相乘
end
function [theta] = exp_ikine(xi_all,g_st,g0)        %反向运动学函数
g_st=g_st/g0;                                       %取消初值部分矩阵
qw=[1.1;0;0.5;1];pw=g_st*qw;pb=[0;0;0.5;1];                         %取3轴上的点作为第一步，化为subp2
theta3=subp3(xi_all(1:4,9:12),qw,pb,norm(pw-pb));
%theta3(1)一脉
[theta1,theta2]=subp2(xi_all(1:4,1:4),xi_all(1:4,5:8),mult(xi_all(1:4,9:12),theta3(1))*qw,pw);
%theta1(1)、theta2(1)一脉
g2=mult([xi_all(1:4,9:12),xi_all(1:4,5:8),xi_all(1:4,1:4)],-[theta3(1) theta2(1) theta1(1)])*g_st;pw2=pb;
[theta4,theta5]=subp2(xi_all(1:4,13:16),xi_all(1:4,17:20),pw2,g2*pw2);pw2=[1;0;0;1];
%theta4(1)、theta5(1)一脉
g3=mult([xi_all(1:4,17:20),xi_all(1:4,13:16),xi_all(1:4,9:12),xi_all(1:4,5:8),xi_all(1:4,1:4)],-[theta5(1) theta4(1) theta3(1) theta2(1) theta1(1)])*g_st;
theta6=subp1(xi_all(1:4,21:24),pw2,g3*pw2);
%theta6(1)
theta=[theta1(1);theta2(1);theta3(1);theta4(1);theta5(1);theta6(1)];%解1
%theta4(2)、theta5(2)一脉
g3=mult([xi_all(1:4,17:20),xi_all(1:4,13:16),xi_all(1:4,9:12),xi_all(1:4,5:8),xi_all(1:4,1:4)],-[theta5(2) theta4(2) theta3(1) theta2(1) theta1(1)])*g_st;
theta6=subp1(xi_all(1:4,21:24),pw2,g3*pw2);
theta=[theta,[theta1(1);theta2(1);theta3(1);theta4(2);theta5(2);theta6(1)]];%解2

%theta1(2)、theta2(2)一脉
g2=mult([xi_all(1:4,9:12),xi_all(1:4,5:8),xi_all(1:4,1:4)],-[theta3(1) theta2(2) theta1(2)])*g_st;pw2=pb;
[theta4,theta5]=subp2(xi_all(1:4,13:16),xi_all(1:4,17:20),pw2,g2*pw2);pw2=[1;0;0;1];
%theta4(1)、theta5(1)一脉
g3=mult([xi_all(1:4,17:20),xi_all(1:4,13:16),xi_all(1:4,9:12),xi_all(1:4,5:8),xi_all(1:4,1:4)],-[theta5(1) theta4(1) theta3(1) theta2(2) theta1(2)])*g_st;
theta6=subp1(xi_all(1:4,21:24),pw2,g3*pw2);
%theta6(1)
theta=[theta,[theta1(2);theta2(2);theta3(1);theta4(1);theta5(1);theta6(1)]];%解3
%theta4(2)、theta5(2)一脉
g3=mult([xi_all(1:4,17:20),xi_all(1:4,13:16),xi_all(1:4,9:12),xi_all(1:4,5:8),xi_all(1:4,1:4)],-[theta5(2) theta4(2) theta3(1) theta2(2) theta1(2)])*g_st;
theta6=subp1(xi_all(1:4,21:24),pw2,g3*pw2);
theta=[theta,[theta1(2);theta2(2);theta3(1);theta4(2);theta5(2);theta6(1)]];%解4

%theta3(2)一脉
theta3(1)=theta3(2);
[theta1,theta2]=subp2(xi_all(1:4,1:4),xi_all(1:4,5:8),mult(xi_all(1:4,9:12),theta3(1))*qw,pw);
%theta1(1)、theta2(1)一脉
g2=mult([xi_all(1:4,9:12),xi_all(1:4,5:8),xi_all(1:4,1:4)],-[theta3(1) theta2(1) theta1(1)])*g_st;pw2=pb;
[theta4,theta5]=subp2(xi_all(1:4,13:16),xi_all(1:4,17:20),pw2,g2*pw2);pw2=[1;0;0;1];
%theta4(1)、theta5(1)一脉
g3=mult([xi_all(1:4,17:20),xi_all(1:4,13:16),xi_all(1:4,9:12),xi_all(1:4,5:8),xi_all(1:4,1:4)],-[theta5(1) theta4(1) theta3(1) theta2(1) theta1(1)])*g_st;
theta6=subp1(xi_all(1:4,21:24),pw2,g3*pw2);
%theta6(1)
theta=[theta,[theta1(1);theta2(1);theta3(1);theta4(1);theta5(1);theta6(1)]];%解5
%theta4(2)、theta5(2)一脉
g3=mult([xi_all(1:4,17:20),xi_all(1:4,13:16),xi_all(1:4,9:12),xi_all(1:4,5:8),xi_all(1:4,1:4)],-[theta5(2) theta4(2) theta3(1) theta2(1) theta1(1)])*g_st;
theta6=subp1(xi_all(1:4,21:24),pw2,g3*pw2);
theta=[theta,[theta1(1);theta2(1);theta3(1);theta4(2);theta5(2);theta6(1)]];%解6

%theta1(2)、theta2(2)一脉
g2=mult([xi_all(1:4,9:12),xi_all(1:4,5:8),xi_all(1:4,1:4)],-[theta3(1) theta2(2) theta1(2)])*g_st;pw2=pb;
[theta4,theta5]=subp2(xi_all(1:4,13:16),xi_all(1:4,17:20),pw2,g2*pw2);pw2=[1;0;0;1];
%theta4(1)、theta5(1)一脉
g3=mult([xi_all(1:4,17:20),xi_all(1:4,13:16),xi_all(1:4,9:12),xi_all(1:4,5:8),xi_all(1:4,1:4)],-[theta5(1) theta4(1) theta3(1) theta2(2) theta1(2)])*g_st;
theta6=subp1(xi_all(1:4,21:24),pw2,g3*pw2);
%theta6(1)
theta=[theta,[theta1(2);theta2(2);theta3(1);theta4(1);theta5(1);theta6(1)]];%解7
%theta4(2)、theta5(2)一脉
g3=mult([xi_all(1:4,17:20),xi_all(1:4,13:16),xi_all(1:4,9:12),xi_all(1:4,5:8),xi_all(1:4,1:4)],-[theta5(2) theta4(2) theta3(1) theta2(2) theta1(2)])*g_st;
theta6=subp1(xi_all(1:4,21:24),pw2,g3*pw2);
theta=[theta,[theta1(2);theta2(2);theta3(1);theta4(2);theta5(2);theta6(1)]];%解8

theta=wrapToPi(theta);         %规整范围到(-pi,pi]
end
function [jacob0,jacobe,gap] = exp_jacob(w_all,q_all,g0,theta_all)
jacob0=zeros(6,size(w_all,2));
xi_all=hight(w_all,q_all);                          %将R^6→se(3)
xi_a=[-cha(w_all,q_all);w_all];
A0=eye(4);
A1=A0*mult(xi_all(1:4,1:4),theta_all(1));
A2=A1*mult(xi_all(1:4,5:8),theta_all(2));
A3=A2*mult(xi_all(1:4,9:12),theta_all(3));
A4=A3*mult(xi_all(1:4,13:16),theta_all(4));
A5=A4*mult(xi_all(1:4,17:20),theta_all(5));
A6=A5*mult(xi_all(1:4,21:24),theta_all(6))*g0;
mat_A=[A0,A1,A2,A3,A4,A5];
for a=1:size(w_all,2)
    jacob0(:,a)=adjoint(mat_A(:,4*(a-1)+1:4*a))*xi_a(:,a);
end
jacobe=adjoint(A6)\jacob0;
gap=gap_mat(A6);
end
function [theta1,theta2] = subp2(xi1,xi2,p,q)        %subp2函数，书上照搬
p=p(1:3);q=q(1:3);
w1=[xi1(3,2);xi1(1,3);xi1(2,1)];w2=[xi2(3,2);xi2(1,3);xi2(2,1)];
r=cal_inter(xi1,xi2);
u=p-r;v=q-r;
alpha=((w1'*w2)*w2'*u-w1'*v)/((w1'*w2)^2-1);
beta=((w1'*w2)*w1'*v-w2'*u)/((w1'*w2)^2-1);
gama=sqrt((norm(u)^2-alpha^2-beta^2-2*alpha*beta*w1'*w2)/norm(xi1(1:3,1:3)*w2)^2);
z1=alpha*w1+beta*w2+gama*(xi1(1:3,1:3)*w2);
z2=alpha*w1+beta*w2-gama*(xi1(1:3,1:3)*w2);
theta1=-[subp1(xi1,q,z1+r),subp1(xi1,q,z2+r)];
theta2=[subp1(xi2,p,z1+r),subp1(xi2,p,z2+r)];
end
function [theta] = subp1(xi,p,q)                  %subp1函数，书上照搬
p=p(1:3);q=q(1:3);
w=[xi(3,2);xi(1,3);xi(2,1)];w_hight=xi(1:3,1:3);v=xi(1:3,4);
r=w_hight*v;
u=p-r;v=q-r;
u1=(eye(3)-w*w')*u;v1=(eye(3)-w*w')*v;
theta=atan2(w'*(cha(u1,v1)),u1'*v1);        
end
function [theta] = subp3(xi,p,q,delta)
p=p(1:3);q=q(1:3);
w=[xi(3,2);xi(1,3);xi(2,1)];w_hight=xi(1:3,1:3);v=xi(1:3,4);
r=w_hight*v;
u=p-r;v=q-r;
u1=(eye(3)-w*w')*u;v1=(eye(3)-w*w')*v;
theta0=atan2(w'*(cha(u1,v1)),u1'*v1);
x=acos((norm(u1)^2+norm(v1)^2-delta^2+abs(w'*(p-q))^2)/(2*norm(u1)*norm(v1)));
theta=[theta0+x,theta0-x];
end
function [adj] = adjoint(A)
a=A(1:3,4);ph=[0 -a(3) a(2);a(3) 0 -a(1);-a(2) a(1) 0];
adj=[A(1:3,1:3),ph*A(1:3,1:3);zeros(3),A(1:3,1:3)];
end
function [gap] = gap_mat(A)
a=A(1:3,4);ph=[0 -a(3) a(2);a(3) 0 -a(1);-a(2) a(1) 0];
gap=[eye(3),-ph;zeros(3),eye(3)];
end
function [r] = cal_inter(xi1,xi2)            %计算两个轴交点的函数，用于解决subp2
w1=[xi1(3,2);xi1(1,3);xi1(2,1)];
gap_p=xi2(1:3,1:3)*xi2(1:3,4)-xi1(1:3,1:3)*xi1(1:3,4);
w_mat=[xi1(3,2) -xi2(3,2);xi1(1,3) -xi2(1,3);xi1(2,1) -xi2(2,1)];
k=w_mat\gap_p;
r=xi1(1:3,1:3)*xi1(1:3,4)+k(1)*w1;
end
function [xi_all] = hight(w_all,q_all)       %%将R^6→se(3)，并补上垂直叉乘部分
xi_all=zeros(4,4*size(w_all,2));
for a=1:size(w_all,2)
    xi_all(1:3,4*(a-1)+1:4*(a-1)+3)=[0,-w_all(3,a),w_all(2,a);w_all(3,a),0,-w_all(1,a);-w_all(2,a),w_all(1,a),0];
    xi_all(1:3,4*a)=-cha(w_all(1:3,a),q_all(1:3,a));
end
end
function [gi_mult] = mult(xi_all,theta_all)
gi_mult=eye(4);
for a=1:length(theta_all)
    w=xi_all(1:3,4*(a-1)+1:4*(a-1)+3);
    v=xi_all(1:3,4*a);
    exp_w=(eye(3)+w*sin(theta_all(a))+w*w*(1-cos(theta_all(a))));
    x=(eye(3)-exp_w)*w*v;
    gi_mult=gi_mult*[exp_w,x;0 0 0 1];
end
end
function [r] = cha(a,b)
r=zeros(3,size(a,2));
for i=1:size(a,2)
    hight_a=[0 -a(3,i) a(2,i);a(3,i) 0 -a(1,i);-a(2,i) a(1,i) 0];
    r(:,i)=hight_a*b(:,i);
end
end