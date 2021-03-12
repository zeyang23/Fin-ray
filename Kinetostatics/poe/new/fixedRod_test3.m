clear
clc

wid=5e-3;
thi=1e-3;
E=197*1e9;

L0=1;
n=50;

pos1=[0.2*L0;0.2*L0;30];

pos2_conv=[0.6*L0;0.6*L0;60];


theta=pos2_conv(3);
end_pos=pos2_conv(1:2);

N=20;

theta_series=linspace(theta,0,N);
endpos_series=[linspace(end_pos(1),L0,N);linspace(end_pos(2),0,N)];


pos2_conv_temp=[endpos_series(:,1);theta_series(1)];    
temp_a=[pos2_conv_temp(1:2);0];
temp_b=rotz(pos1(3))*temp_a;
pos2=[temp_b(1:2)+pos1(1:2);pos2_conv_temp(3)+pos1(3)];
Rod= fixedRod(E,L0,wid,thi,n,pos1,pos2);
Rod.init_exp;


for i=1:N
    theta_last=Rod.conv_theta;
    
    pos2_conv_temp=[endpos_series(:,i);theta_series(i)];
    
    temp_a=[pos2_conv_temp(1:2);0];
    temp_b=rotz(pos1(3))*temp_a;

    pos2=[temp_b(1:2)+pos1(1:2);pos2_conv_temp(3)+pos1(3)];
    
    Rod = fixedRod(E,L0,wid,thi,n,pos1,pos2);
    Rod.init_exp;
    Rod.conv_theta=theta_last;
    Rod.update_conv;
    
    TOL=1e-6;
    Rod.Newton_conv(TOL);

    Rod.plot_pos;
end