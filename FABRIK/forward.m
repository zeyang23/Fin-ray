% 总共有n个点，target是终点坐标 1*2
% pi是n个点上一步的坐标 n*2
% d是n个点之间的距离 (n-1)*1
% pi_new是新返回的点的坐标 n*2
function p_new=forward(p,d,target)
    n=length(d)+1;
    p_new=p;
    p_new(end,:)=target;
    
    for i=n-1:-1:1
        ri=norm(p_new(i+1,:)-p_new(i,:));
        di=d(i);
        lambda=di/ri;
        p_new(i,:)=(1-lambda)*p_new(i+1,:)+lambda*p_new(i,:);
    end
end