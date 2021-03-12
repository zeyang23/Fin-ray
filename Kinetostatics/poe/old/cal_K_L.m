function K_L=cal_K_L(I,E,n,x)
    L=x(1);
    F=x(end-5:end);
    theta=x(2:end-6);

    %΢Ԫ�ĳ���
    delta=L/n;

    %���ϵ��������
    c=delta/(E*I);

    %�նȾ���
    K=diag((1/c)*ones(n,1));
    
    K_L=-K*theta'/L;
end