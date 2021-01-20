function draw_finray(theta1,theta2,theta3)
    L1=400;
    L2=200/cos(14/180*pi);
    phi=14/180*pi;
    psi=2*(pi/2-phi);

    L3=500;

    A0=[700,0];

    A1=[A0(1)-L1*sin(theta1),A0(2)+L1*cos(theta1)];
    A2=[A1(1)-L1*sin(theta1+theta2),A1(2)+L1*cos(theta1+theta2)];
    A3=[A2(1)-L1*sin(theta1+theta2+theta3),A2(2)+L1*cos(theta1+theta2+theta3)];
    draw_line(A0,A1)
    hold on 
    axis equal
    axis([0 800 0 1200])
    draw_line(A1,A2)
    draw_line(A2,A3)

    B1=[A0(1)-L2*sin(theta1+phi),A0(2)+L2*cos(theta1+phi)];
    draw_line(A0,B1)
    draw_line(B1,A1)

    B2=[A1(1)-L2*sin(theta1+theta2+phi),A1(2)+L2*cos(theta1+theta2+phi)];
    draw_line(A1,B2)
    draw_line(B2,A2)

    B3=[A2(1)-L2*sin(theta1+theta2+theta3+phi),A2(2)+L2*cos(theta1+theta2+theta3+phi)];
    draw_line(A2,B3)
    draw_line(B3,A3)

    P0=[0,0];
    C1=solve_triangle(P0,B1,L2,L3);
    draw_line(P0,C1)
    draw_line(C1,B1)

    P1=solve_triangle2(P0,C1,psi);
    draw_line(C1,P1)
    draw_line(P0,P1)

    C2=solve_triangle(P1,B2,L2,300);
    draw_line(P1,C2)

    P2=solve_triangle2(P1,C2,psi);
    draw_line(C2,P2)
    draw_line(P1,P2)

    draw_line(C2,B2)

    C3=solve_triangle(P2,B3,L2,100);
    draw_line(P2,C3)

    P3=solve_triangle2(P2,C3,psi);
    draw_line(C3,P3)
    draw_line(P3,P3)

    draw_line(P2,P3)
    draw_line(C3,B3)

    draw_line(P3,A3)
end