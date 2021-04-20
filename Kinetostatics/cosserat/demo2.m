% 末端受力，求解形状
% shooting method

clear
clc

wid=14e-3;
thi=0.3e-3;
E=197*1e9;

Iz=1/12*thi.^3*wid;

L=1;

Fe=[-0.048229138674937;-0.048229143826622;-6.952354689790895e-04];

f=@(x) check_balance(x,L,E,Iz,Fe);

x0=zeros(3,1);

options2 = optimoptions('fsolve','Display','off');
x_cosserat=fsolve(f,x0,options2);


% 验证结果

span = [0 L];
y0 = [0;0;0;x_cosserat];
options3=odeset('MaxStep',1e-2);
[s,Y] = ode45(@(s,y) get_ydot(s,y,L,E,Iz,Fe), span, y0,options3);

hold on
plot(Y(:,1),Y(:,2),'green')



function res=check_balance(x,L,E,I,Fe)
    
    n0=x(1:2);
    m0=x(3);

    span = [0 L];
    y0 = [0;0;0;n0;m0];
    
    options=odeset('MaxStep',1e-2);
    
    [~,Y] = ode45(@(s,y) get_ydot(s,y,L,E,I,Fe), span, y0,options);
    
    ye=transpose(Y(end,:));
    
    res=zeros(3,1);
    res(1:2)=ye(4:5)-Fe(1:2);
    res(3)=transpose(ye(1:2))*[0 1;-1 0]*ye(4:5)+ye(6)-transpose(ye(1:2))*[0 1;-1 0]*Fe(1:2)-Fe(3);
end


function ydot=get_ydot(s,y,L,E,I,Fe)
    ydot=zeros(size(y));
    
    r=y(1:2);
    theta=y(3);
    n=y(4:5);
    m=y(6);
    
%     delta=1e-16;
%     if (L-s)<delta
%         f=Fe(1:2);
%         l=Fe(3);
%     else
%         f=[0;0];
%         l=0;
%     end

    f=[0;0];
    l=0;
    
    ydot(1:2)=[cos(theta);sin(theta)];
    ydot(3)=1/(E*I)*m;
    ydot(4:5)=-f;
    ydot(6)=-[cos(theta) sin(theta)]*[0 1;-1 0]*n-l;
end