function JACOBIANs=Jacob_constraint(I,E,n,x)
    % x: L F theta
    L=x(1);
    F=x(end-5:end);
    theta=x(2:end-6);

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
    
    jacobs=exp_jacob_simple(w_all,q_all,theta);
    
    
    J_L=cal_J_L(w_all,q_all,x);
    K_L=cal_K_L(I,E,n,x);
    K_J=-partial_J_theta(w_all,q_all,theta,jacobs,F');
    
    JACOBIANs=zeros(length(x));
    JACOBIANs(1:6,1)=J_L;
    JACOBIANs(1:6,2:end-6)=jacobs;
    JACOBIANs(7:end-1,1)=K_L;
    JACOBIANs(7:end-1,2:end-6)=K+K_J;
    JACOBIANs(7:end-1,end-5:end)=-transpose(jacobs);
    JACOBIANs(end,end-5:end)=[1 0 0 0 0 0];
end