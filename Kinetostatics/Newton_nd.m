function xsolve=Newton_nd(f,Jacob,x0,TOL)
    % f返回n*1列向量 
    % Jacob返回n*n矩阵
    % x是n*1列向量
    
    x(:,1)=x0';
    k=1;
    while(1)
        if k>500
            error("can not converge")
        end
        J=Jacob(x(:,k)');
        b=f(x(:,k)');
        delta=-pinv(J)*b;
        x(:,k+1)=x(:,k)+delta;
        k=k+1;
%         if(norm(x(:,k)-x(:,k-1))<TOL)
        if(norm(f(x(:,k)'))<TOL)
            fprintf('Newton Method converge: iteration = %d\n',k-1)
            fprintf('norm(e) = %E\n',norm(f(x(:,k)')))
            break;
        end
    end
    xsolve=x(:,end);
    xsolve=xsolve';
end