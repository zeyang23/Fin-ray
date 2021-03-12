% 21-03-08 ����ʹ�ö��������ķ�����⻬��Լ�������
% ������ѵ㣺�������ɺ��ʵ����ֵ����һ��������������Ҫ���ַ����������ʵĳ�ֵ

% �����ƣ�
% 1 ��дǰ����ģ�飬���ɺ��ʵĳ�ֵ
% 2 �����Ƶ�֤������㸽����f(l)�ǵ�����

clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;

pos1=[0*L0;0*L0;0];

pos2=[0.3*L0;0.3*L0;60];

dist=sqrt((pos1(1)-pos2(1))^2+(pos1(2)-pos2(2))^2);

La=1.05*dist;
Lb=1.4*dist;

TOL=0.5e-4; %��ȷ��С�����6λ

f=@(L) func(L,pos2);

Lsolve=Bisection(f,La,Lb,TOL);


RodA=fixedRod(E,Lsolve,wid,thi,n,pos1,pos2);
RodA.init_exp;
TOL=1e-6;
RodA.Newton_conv(TOL);
RodA.plot_pos;

pos3=[0.6*L0;0.6*L0;90];

RodB=fixedRod(E,L0-Lsolve,wid,thi,n,pos2,pos3);
RodB.init_exp;
TOL=1e-6;
RodB.Newton_conv(TOL);
RodB.plot_pos;


% �������е�����Բ
pos2theta=pos2(3);
pos2x=pos2(1);
pos2y=pos2(2);

r=0.01*L0;
vec=r*[cos(pos2theta/180*pi);sin(pos2theta/180*pi);0];
vec1=rotz(90)*vec;
vec2=rotz(-90)*vec;
center1x=vec1(1)+pos2x;
center1y=vec1(2)+pos2y;
center2x=vec2(1)+pos2x;
center2y=vec2(2)+pos2y;

rectangle('Position',[center1x-r,center1y-r,2*r,2*r],'Curvature',[1,1],'linewidth',1,'edgecolor','r')
rectangle('Position',[center2x-r,center2y-r,2*r,2*r],'Curvature',[1,1],'linewidth',1,'edgecolor','r')

function r=func(L,pos2)
    wid=5e-3;
    thi=1e-3;
    E=197*1e9;
    L0=1;

    pos1=[0*L0;0*L0;0];
    
    Rod = fixedRod(E,L,wid,thi,50,pos1,pos2);

    Rod.init_exp;
    TOL=1e-6;
    Rod.Newton_conv(TOL);
    
    theta=pos2(3)/180*pi;
    dir=[cos(theta);sin(theta);0;0;0;0];

    r=transpose(dir)*Rod.conv_F;

end

function x=Bisection(f,a,b,TOL) 
    % f�Ǻ���������������������
    % a��������ˣ�b�������Ҷ�
    % TOL=0.5e-p ��ȷ��С�����pλ
    
    % ����a��b�Ƿ���������
    if(sign(f(a))==sign(f(b)))
        error('wrong input');
    end
    
    while (b-a)/2 > TOL
        c=(a+b)/2;
        fc=f(c);
        if(fc==0)
            break;
        end
        if sign(fc)==sign(f(b))
            b=c;
        else
            a=c;
        end
    end
    x=(a+b)/2;
end