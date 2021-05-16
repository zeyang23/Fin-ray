function Blist=SlistToBlist(Slist,M)
    Blist=zeros(size(Slist));
    for i=1:size(Blist,2)
        Blist(:,i)=Adjoint(TransInv(M))*Slist(:,i);
    end
end