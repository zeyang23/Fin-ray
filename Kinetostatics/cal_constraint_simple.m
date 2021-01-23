function r=cal_constraint_simple(L,I,E,n,gt,x)
    % x: theta F
    F=x(end-5:end);
    theta=x(1:end-6);

    %微元的长度
    delta=L/n;

    %材料的物理参数
    c=delta/(E*I);

    %刚度矩阵
    K=diag((1/c)*ones(n,1));
    
    %所有的w向量构成矩阵
    w_all=[];
    for i=1:n
        w_all=[w_all,[0;0;1]];
    end

    %所有的q向量构成矩阵
    q_all=zeros(3,n);
    for i=1:n
        q_all(:,i)=[1;0;0]*delta*(i-1/2);
    end

    %初始姿态矩阵
    g0=[eye(3),[L;0;0];0,0,0,1];
    
    jacobs=exp_jacob_simple(w_all,q_all,theta);
    [g,~]=exp_fkine(w_all,q_all,g0,theta); 
    
    tau=K*theta'-transpose(jacobs)*F';
    
    error=logm(inv(gt)*g);
    e(1:3,1)=error(1:3,4);
    e(4,1)=-error(2,3);
    e(5,1)=error(1,3);
    e(6,1)=-error(1,2);
    
    fx=[1;0;0;0;0;0]'*F';
    
    r=[e;tau];
    
end