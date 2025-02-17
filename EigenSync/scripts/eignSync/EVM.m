function rots = EVM (M,E,w)


n= length(M)/3;


G=sparse([E(:,1)',E(:,2)'],[E(:,2)',E(:,1)'],[w,w],n,n);
D=sparse(1:3*n,1:3*n,1,3*n,3*n);
d=sum(G);

for i= 1:n
    D(3*i-2:3*i,3*i-2:3*i)=eye(3)*d(i);
end
D=sparse(D);
M_norm=D^(-1)*M;

[v,e]=eigs(M_norm,3,'lr');

v=v(:,1:3);

rots=cell(n,1);


for i=1:n
    x=v((3*i-2):(3*i),1:3);
    [U,~,V]=svd(x);
    rots{i}=V*U';
end

end