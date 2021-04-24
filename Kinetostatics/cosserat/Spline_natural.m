function v=Spline_natural(x,y,u)
    delta_x=diff(x);
    delta_y=diff(y);
    n=length(x);
    A=zeros(n);
    A(1)=1;
    A(end)=1;
    for ii=2:n-1
        A(ii,ii-1)=delta_x(ii-1);
        A(ii,ii)=2*(delta_x(ii-1)+delta_x(ii));
        A(ii,ii+1)=delta_x(ii);
    end
    B=zeros(n,1);
    B(1)=0;
    B(n)=0;
    for jj=2:n-1
        B(jj)=3*(delta_y(jj)./delta_x(jj)-delta_y(jj-1)./delta_x(jj-1));
    end
    
    % 生成系数 a b c d 均为1×(n-1)行向量
    c=A\B;
    c=transpose(c);
    b=zeros(1,n-1);
    for kk=1:length(b)
        b(kk)=delta_y(kk)./delta_x(kk)-2/3*delta_x(kk)*c(kk)-1/3*delta_x(kk)*c(kk+1);
    end
    d=zeros(1,n-1);
    for kk=1:length(d)
        d(kk)=1./(3*delta_x(kk))*(c(kk+1)-c(kk));
    end
    c(end)=[];
    a=zeros(1,n-1);
    for kk=1:length(a)
        a(kk)=y(kk);
    end
    
    % 根据u计算v 利用了数组的逻辑索引法
    k=ones(size(u));
    for ii=2:n-1
        k(x(ii)<=u)=ii;
    end
    s=u-x(k);
    v=a(k)+s.*(b(k)+s.*(c(k)+s.*d(k)));
end
    