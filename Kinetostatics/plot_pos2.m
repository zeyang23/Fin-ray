function plot_pos2(w_all,q_all,g0,delta,theta_solve,L0)
    theta=zeros(1,length(theta_solve));
    for k=1:length(theta)
        theta(k)=theta_solve(k);
    end
    
    [g_exp,~]=exp_fkine(w_all,q_all,g0,theta);

    %得到整个杆的位姿
    pos=[0;0];
    pos=[pos,[q_all(1,1);q_all(2,1)]];
    for i =2:length(theta)
        g0i=[eye(3),[(i-0.5)*delta;0;0];0,0,0,1];
        [gsti,~]=exp_fkine(w_all(:,1:i-1),q_all(:,1:i-1),g0i,theta(:,1:i-1));
        pos=[pos,gsti(1:2,4)];
    end
    pos=[pos,g_exp(1:2,4)];

    plot(pos(1,:),pos(2,:),'-o','MarkerSize',2)
    
    
    new_xy=g_exp(1:2,1:2)*[1 0;0 1];
    new_x=new_xy(1,:);
    new_y=new_xy(2,:);
    hold on
    end_pos=g_exp(1:2,4);
    location_x=[end_pos(1),end_pos(1)];
    location_y=[end_pos(2),end_pos(2)];
    quiver(location_x,location_y,new_x,new_y,0.05)
    axis equal
    
    axis([0,L0,-0.2*L0,0.6*L0])
end