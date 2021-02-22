% 总共有n个点，base是起点坐标 1*2
% pi是n个点上一步的坐标 n*2
% d是n个点之间的距离 (n-1)*1
% pi_new是新返回的点的坐标 n*2
function p_new=backward(p,d,base)
    n=length(d)+1;
    p_new=p;
    p_new(1,:)=base;
    
    for i=1:n-1
        ri=norm(p_new(i+1,:)-p_new(i,:));
        di=d(i);
        lambda=di/ri;
        p_new(i+1,:)=(1-lambda)*p_new(i,:)+lambda*p_new(i+1,:);
    end
end