function JACOBIANs=Jacob_constraint_simple2(L,I,E,n,x)
    % x: theta F
    F=x(end-5:end);
    theta=x(1:end-6);

    %΢Ԫ�ĳ���
    delta=L/n;

    %���ϵ��������
    c=delta/(E*I);

    %�նȾ���
    K=diag((1/c)*ones(n,1));
    
    %���е�w�������ɾ���
    w_all=[];
    for i=1:n
        w_all=[w_all,[0;0;1]];
    end


    %���е�q�������ɾ���
    q_all=zeros(3,n);
    for i=1:n
        q_all(:,i)=[1;0;0]*delta*(i-1/2);
    end

    %��ʼ��̬����
    g0=[eye(3),[L;0;0];0,0,0,1];
    
    [jacobs,jacobe,~]=exp_jacob(w_all,q_all,g0,theta);
    
    K_J=-partial_J_theta(w_all,q_all,theta,jacobs,F');
    
    JACOBIANs=zeros(length(x));
    JACOBIANs(1:6,1:n)=-jacobe;
    JACOBIANs(7:end,1:n)=K+K_J;
    JACOBIANs(7:end,n+1:end)=-transpose(jacobs);
end