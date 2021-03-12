function tau=cal_balance(K,w_all,q_all,g0,theta_all,F)
    [jacobs,~,~]=exp_jacob(w_all,q_all,g0,theta_all);
    tau=K*transpose(theta_all)-transpose(jacobs)*F;
end