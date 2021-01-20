clear
clc

theta1=-6/180*pi;
theta2=15/180*pi;
theta3=30/180*pi;

% draw_finray(14.5/180*pi,0,0)

n=100;

theta1_serie=linspace(14.5/180*pi,theta1,n);
theta2_serie=linspace(0,theta2,n);
theta3_serie=linspace(0,theta3,n);

for i=1:n
    thetaA=theta1_serie(i);
    thetaB=theta2_serie(i);
    thetaC=theta3_serie(i);
    
    draw_finray(thetaA,thetaB,thetaC)
    pause(0.01)
    if i<n
        clf
    end
end