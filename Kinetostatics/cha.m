function [r] = cha(a,b)
    r=zeros(3,size(a,2));
    for i=1:size(a,2)
        hight_a=[0 -a(3,i) a(2,i);a(3,i) 0 -a(1,i);-a(2,i) a(1,i) 0];
        r(:,i)=hight_a*b(:,i);
    end
end