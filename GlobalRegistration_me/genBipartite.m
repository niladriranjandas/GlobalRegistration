n=20; %number of points
m=5;  %number of patches
d=3;  %dimenssion

edges=[];
nOutDeg = zeros(1,n);
for i=1:n
    nOutDeg(i) = randi(m);
    
    endpts = sort(randperm(m,nOutDeg(i)));
    edges_ = [repmat(i,nOutDeg(i),1),endpts'];
    edges = [edges;edges_];
end
    
clearvars edges_    


%%%%%%%%%%%%%%%%%%%%
x_k_i = randi([1,50],3,size(edges,1));

%%%% B %%%%%%
B = zeros(m*d,n+m);
tmpB = zeros(d,n+m);
i=2; % i denotes patch-i

for i=1:m
    
     e_m_i = zeros(m,1);
     e_m_i(i) = 1;     
        
     e_kron_id = kron(e_m_i,eye(d));
     index = find(edges(:,2)==i);
   
     x_k_i_ = x_k_i(:,index);
     pts = edges(index,1);
     patch = i;
   
     for j=1:size(x_k_i_,2)
         e_k_i  = zeros(n+m,1);
         e_k_i(pts(j)) = 1;e_k_i(n+i) = -1;

         B = B + e_kron_id * x_k_i_(:,j) * e_k_i';
        tmpB = tmpB + x_k_i_(:,j) * e_k_i';
     end
end    

%%%%%%%%% adjacenecy matrix and L %%%%%%%

% for i=1:size(edges,1)
%         adjacent(edges(i,1),n+edges(i,2)) = 1;
%         adjacent(n+edges(i,2),edges(i,1)) = 1;
% end
  adj_mat = sparse(edges(:,1),n+edges(:,2),1,m+n,m+n);
  adj_mat = adj_mat + adj_mat';
  
  deg_mat = sparse(1:m+n,1:m+n,sum(adj_mat),m+n,m+n);
  L_calc = deg_mat - adj_mat;
  
%   L = sparse(n+m,n+m);
%   tmpL = zeros(n+m);
%   for i=1:m
%      index = find(edges(:,2)==i);
%      pts = edges(index,1);
%       for j = 1:size(index,1)
%           e_k_i = zeros(m+n);
%           e_k_i(pts(j)) = 1; e_k_i(n+i) = -1;
%           
%           L = L+ e_k_i * e_k_i';
%           tmpL = tmpL + e_k_i * e_k_i';
%       end
%       tmpL
%   end

%%%%%% D %%%%%%%%%
D = zeros(m*d);
tmpD = zeros(d);

i=3;
   
   index = find(edges(:,2)==i);
   x_k_i_ = x_k_i(:,index);
   patch=i;
   
   e_m_i = zeros(m,1);
   e_m_i(i) = 1;
   for j=1:size(x_k_i_,2)
       e_kron_id = kron(e_m_i,eye(d));
       D = D + e_kron_id * x_k_i_(:,j) * x_k_i_(:,j)' * e_kron_id';
       tmpD = tmpD + x_k_i_(:,j) * x_k_i_(:,j)';
   end
   