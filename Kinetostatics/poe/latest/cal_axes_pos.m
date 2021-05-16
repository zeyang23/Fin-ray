function pos_all=cal_axes_pos(thetalist,Qlist,Slist)
    n=length(thetalist);
    pos_all=zeros(2,n);
    
    for i=1:n
        M=[eye(3),Qlist(:,i);[0 0 0 1]];
        T=FKinSpace(M, Slist(:,1:i), thetalist(1:i));
        
        [~,p]=TransToRp(T);
        
        pos_all(:,i)=[p(1);p(2)];
    end
end