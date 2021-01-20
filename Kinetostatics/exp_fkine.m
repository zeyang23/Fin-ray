function [g_st,xi_all] = exp_fkine(w_all,q_all,g0,theta_all)%正向运动学函数
    xi_all=hight(w_all,q_all);                          %将R^6→se(3)
    gi_mult=mult(xi_all,theta_all);                     %将se(3)→SE(3)并累乘
    g_st=gi_mult*g0;                                    %初值矩阵相乘
end