function theta=p2theta(p)
    % p��n*2���� ����n+2�����λ�� ��� n���ؽ�λ�� �յ�
    % theta n*1������
    
    % �����õ�n+1���߶ε�x y����
    X=diff(p(:,1));
    Y=diff(p(:,2));
    
    n=length(X)-1;
    
    theta=zeros(n,1);
    
    for i=1:n
        thetaA=atan2(Y(i),X(i));
        thetaB=atan2(Y(i+1),X(i+1));
        theta(i)=thetaB-thetaA;
    end
end