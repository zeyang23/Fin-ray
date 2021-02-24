function theta=p2theta(p)
    % p是n*2向量 代表n+2个点的位置 起点 n个关节位置 终点
    % theta n*1列向量
    
    % 首先拿到n+1个线段的x y坐标
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