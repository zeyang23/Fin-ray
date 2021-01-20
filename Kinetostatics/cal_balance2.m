function tau=cal_balance(K,w_all,q_all,g0,theta_all,F)
    [jacobs,~,~]=exp_jacob(w_all,q_all,g0,theta_all);
    j_up=jacobs(1:3,:);
    j_down=jacobs(4:6,:);
    J=[j_down;j_up];
    tau=K*transpose(theta_all)-transpose(J)*F;
end