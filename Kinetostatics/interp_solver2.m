function x=interp_solver2(theta,end_pos,N,n,L0,I,E,Jacob)
    
    theta_series=linspace(0,theta,N);
    endpos_series=[linspace(L0,end_pos(1),N);linspace(0,end_pos(2),N)];
    
    x0=zeros(n+6,1);
    TOL=1e-6;
    
    for i=1:N
        THETA=theta_series(i);
        ENDPOS=endpos_series(:,i);
        gt=[rotz(THETA/180*pi),[ENDPOS;0];0,0,0,1];
        
        f=@(x) cal_constraint_simple2(L0,I,E,n,gt,x);
        
        xsolve=Newton_nd(f,Jacob,x0,TOL);
        x0=xsolve;
    end
    
    x=x0;
end