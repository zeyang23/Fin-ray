function C=solve_triangle(A,B,L1,L2)
    AB=B-A;
    theta=atan2(AB(2),AB(1));
    L3=norm(B-A);
    thetaA=acos((L1^2+L3^2-L2^2)/(2*L1*L3));
    thetaC=theta+thetaA;
    C=A+[L1*cos(thetaC),L1*sin(thetaC)];
end
    
    