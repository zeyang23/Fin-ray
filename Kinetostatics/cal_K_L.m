function K_L=cal_K_L(I,E,n,x)
    L=x(1);
    F=x(end-5:end);
    theta=x(2:end-6);

    %微元的长度
    delta=L/n;

    %材料的物理参数
    c=delta/(E*I);

    %刚度矩阵
    K=diag((1/c)*ones(n,1));
    
    K_L=-K*theta'/L;
end