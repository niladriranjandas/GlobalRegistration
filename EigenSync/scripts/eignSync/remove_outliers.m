function E_new = remove_outliers(E,d)


cutoff=mean(d)+std(d);

O=find(d>cutoff);

n=max(E(:,2));
E_new = E;

E_new(O,:)=[];
G=sparse(E_new(:,2),E_new(:,1),1,n,n);
[num_comp,C]=graphconncomp(G,'Directed','false');

if(num_comp==1)
    return
end


E_out=C(E(O,:))';d_out=d(O);D=O;

same_comp=find(E_out(:,1)-E_out(:,2)==0);
E_out(same_comp,:)=[];
d_out(same_comp)=[];D(same_comp)=[];

x=sort(E_out,2);x=sortrows([x,d_out,D]);
[b,m,~]=unique(x(:,1:2),'rows','first');disp(num_comp);disp(b);
comp_graph=sparse(b(:,1),b(:,2),x(m,3),num_comp,num_comp);
MST=graphminspantree(comp_graph);

[x,y,~]=find(MST);
[~,~,idx]=intersect([y,x],[b(:,1),b(:,2)],'rows');
O_new=setdiff(O,D(idx));
E_new=E;
E_new(O_new,:)=[];


end