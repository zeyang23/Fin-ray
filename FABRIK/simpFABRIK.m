function p_solve=simpFABRIK(base,target,d,p0,tol)
    k=0;
    
    p=p0;
    while(1)
        if norm(p(end,:)-target)<tol
            fprintf('FABRIK converge: iteration = %d\n',k)
            break;
        elseif k>1000
            error('fail to converge')
        else
            p=forward(p,d,target);
            p=backward(p,d,base);

            k=k+1;
        end
    end
    p_solve=p;
end