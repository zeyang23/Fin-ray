% 21-05-16

% ���ݶȴ���fsolve����Ȼ�޷����
% ���ң�����ݶ��޷�ͨ��fsolve���ݶȼ��
% ��ˣ��Ҿ��û���poe�����˶�ѧ�㷨����ݶ��ǽ���ֵ


clear
clc


wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;

pos1=[0*L0;0*L0;0];
pos2=[0.6*L0;0.6*L0;60];
RodA = fixedRod(E,L0,wid,thi,n,pos1,pos2);

RodA.init_exp;
% TOL=1e-6;
% RodA.Newton_conv(TOL);

x0=zeros(n+6,1);
f=@(x) cal_balance(x,RodA);

options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false);
[x_solve,fval,exitflag,output] = fsolve(f,x0,options);


RodA.plot_pos_conv;

function [r,J]=cal_balance(x,Rod)
    Rod.conv_theta=x(1:Rod.n_seg);
    Rod.conv_F=x(Rod.n_seg+1:Rod.n_seg+6);
    
    Rod.update_conv;
    
    r=Rod.cal_constraint_conv;
    J=Rod.cal_Jacob_conv;
                
end