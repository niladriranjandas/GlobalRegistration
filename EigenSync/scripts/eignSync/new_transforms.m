function [new_M,new_T,new_w]=new_transforms(M,T,E,w,E_sub)



[~,~,ids]=intersect(E_sub,E,'rows');  %disp(size(T));disp(ids);
new_T=T(ids,:);
new_w=w(ids);
new_M=spalloc(length(M),length(M),18*length(E_sub));
for k=ids
    i=E(k,1);
    j=E(k,2);
    new_M(3*i-2:3*i,3*j-2:3*j)= M(3*i-2:3*i,3*j-2:3*j);
    new_M(3*i-2:3*i,3*j-2:3*j)= M(3*i-2:3*i,3*j-2:3*j);
end









end