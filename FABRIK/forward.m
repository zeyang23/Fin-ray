% �ܹ���n���㣬target���յ����� 1*2
% pi��n������һ�������� n*2
% d��n����֮��ľ��� (n-1)*1
% pi_new���·��صĵ������ n*2
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