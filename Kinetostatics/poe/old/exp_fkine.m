function [g_st,xi_all] = exp_fkine(w_all,q_all,g0,theta_all)%�����˶�ѧ����
    xi_all=hight(w_all,q_all);                          %��R^6��se(3)
    gi_mult=mult(xi_all,theta_all);                     %��se(3)��SE(3)���۳�
    g_st=gi_mult*g0;                                    %��ֵ�������
end