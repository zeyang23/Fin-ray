function plot_pos(w_all,q_all,g0,delta,theta_solve)
    theta=zeros(1,length(theta_solve));
    for k=1:length(theta)
        theta(k)=theta_solve(k);
    end
    
    [g_exp,~]=exp_fkine(w_all,q_all,g0,theta);

    %得到整个杆的位姿
    pos=[0;0];
    pos=[pos,[q_all(1,2);q_all(2,2)]];
    for i =2:length(theta)
        g0i=[eye(3),[(i-0.5)*delta;0;0];0,0,0,1];
        [gsti,~]=exp_fkine(w_all(:,1:i-1),q_all(:,1:i-1),g0i,theta(:,1:i-1));
        pos=[pos,gsti(1:2,4)];
    end
    pos=[pos,g_exp(1:2,4)];

    plot(pos(1,:),pos(2,:),'-o','MarkerSize',2)
end