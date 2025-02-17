function [rots,trans]   =   TSR_global_align(M,T,E,w,num_iter)


E1=E;
M1=M;
T1=T;
w1=w;

for n = 1:num_iter
    rots=EVM(M1,E1,w1);
    m=length(E1);
    d=zeros(E,1);
    for k=1:m
        i = E1(k,1); j = E1(k,2);
        actual = full(M1(3*i-2:3*i,3*j-2:3*j));
        estim =rots{i}*rots{j}';
        d(k)= norm(actual-estim,'fro');
    end
    E1=remove_outliers(E1,d);
    [M1,T1,w1]=new_transforms(M,T,E,w,E1);
end

rots = EVM(M1,E1,w1);
 trans = recover_trans(T1,E1,rots);







end