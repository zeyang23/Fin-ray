clear
clc

load('A_start.mat')
load('B_start.mat')
load('C_start.mat')
load('D_start.mat')

load('A_end.mat')
load('B_end.mat')
load('C_end.mat')
load('D_end.mat')

N=21;

for k=1:52
    Ax_series(k,:)=linspace(A_start(1,k),A_end(1,k),N);
    Ay_series(k,:)=linspace(A_start(2,k),A_end(2,k),N);

    Bx_series(k,:)=linspace(B_start(1,k),B_end(1,k),N);
    By_series(k,:)=linspace(B_start(2,k),B_end(2,k),N);

    Cx_series(k,:)=linspace(C_start(1,k),D_end(1,k),N);
    Cy_series(k,:)=linspace(C_start(2,k),D_end(2,k),N);

    Dx_series(k,:)=linspace(D_start(1,k),C_end(1,k),N);
    Dy_series(k,:)=linspace(D_start(2,k),C_end(2,k),N);
end

avi_object = VideoWriter('fake.avi');
avi_object.FrameRate = 20;
open(avi_object);

figure
for I=1:N
    hold on
    axis equal
    
    Ax=Ax_series(:,I);
    Ay=Ay_series(:,I);
    plot(Ax,Ay,'b-o','MarkerSize',2)
    
    Bx=Bx_series(:,I);
    By=By_series(:,I);
    plot(Bx,By,'b-o','MarkerSize',2)
    
    Cx=Cx_series(:,I);
    Cy=Cy_series(:,I);
    plot(Cx,Cy,'b-o','MarkerSize',2)
    
    Dx=Dx_series(:,I);
    Dy=Dy_series(:,I);
    plot(Dx,Dy,'b-o','MarkerSize',2)
    
    A0=[Ax_series(1,I),Ay_series(1,I)];
    A1=[Ax_series(8,I),Ay_series(8,I)];
    A2=[Ax_series(16,I),Ay_series(16,I)];
    A3=[Ax_series(24,I),Ay_series(24,I)];
    A4=[Ax_series(32,I),Ay_series(32,I)];
    A5=[Ax_series(40,I),Ay_series(40,I)];
    
    B0=[Bx_series(1,I),By_series(1,I)];
    B1=[Bx_series(8,I),By_series(8,I)];
    B2=[Bx_series(16,I),By_series(16,I)];
    B3=[Bx_series(24,I),By_series(24,I)];
    B4=[Bx_series(32,I),By_series(32,I)];
    B5=[Bx_series(40,I),By_series(40,I)];
    
    C0=[Cx_series(1,I),Cy_series(1,I)];
    C1=[Cx_series(8,I),Cy_series(8,I)];
    C2=[Cx_series(16,I),Cy_series(16,I)];
    C3=[Cx_series(24,I),Cy_series(24,I)];
    C4=[Cx_series(32,I),Cy_series(32,I)];
    C5=[Cx_series(40,I),Cy_series(40,I)];
        
    D0=[Dx_series(1,I),Dy_series(1,I)];
    D1=[Dx_series(8,I),Dy_series(8,I)];
    D2=[Dx_series(16,I),Dy_series(16,I)];
    D3=[Dx_series(24,I),Dy_series(24,I)];
    D4=[Dx_series(32,I),Dy_series(32,I)];
    D5=[Dx_series(40,I),Dy_series(40,I)];
    
    draw_line(C0,D0)
    draw_line(C1,D1)
    draw_line(C2,D2)
    draw_line(C3,D3)
    draw_line(C4,D4)
    draw_line(C5,D5)

    draw_line(A0,B0)
    draw_line(A1,B1)
    draw_line(A2,B2)
    draw_line(A3,B3)
    draw_line(A4,B4)
    draw_line(A5,B5)
    
    object_right=B_end(:,20:42);
    object_left=D_end(:,20:42);

    for i =1:22
        A=object_right(:,i);
        A(1)=A(1)-0.01;
        B=object_right(:,i+1);
        B(1)=B(1)-0.01;
        C=object_left(:,i);
        C(1)=C(1)+0.01;
        D=object_left(:,i+1);
        D(1)=D(1)+0.01;
        draw_line2(A,B)
        draw_line2(C,D)
    end

    R1=object_right(:,end);
    R1(1)=R1(1)-0.01;
    R2=object_right(:,1);
    R2(1)=R2(1)-0.01;
    L1=object_left(:,end);
    L1(1)=L1(1)+0.01;
    L2=object_left(:,1);
    L2(1)=L2(1)+0.01;
    draw_line2(L1,R1)
    draw_line2(L2,R2)
    
    axis([-1 0.4 -0.1 1])

    M = getframe;
    writeVideo(avi_object,M);

    if I < N
        clf;
    end
end
close(avi_object);