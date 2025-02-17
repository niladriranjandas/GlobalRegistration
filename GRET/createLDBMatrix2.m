function [L,D,B]=createLDBMatrix2(patches,N)


disp('Creating L D B matix');
M   = size(patches,1);
dim = size(patches{1,2},1);
% set up the block matrices

%L = sparse(zeros(N + M, N + M));
B = sparse(zeros(M*dim, N + M));
D = sparse(zeros(M*dim, M*dim));;
%% L matrix
% structure of L matrix 
% L= (temp_1)_(N x N) | (temp_2)_(N x M)
%    --------------------
%    (temp_3)_(M x N) | (temp_4)_(M x M)
% 
% temp_1 = diagonal matrix  number of patches each point is 
% temp_3 and temp_2 = position of point in each patch
% temp_4 = diagonal of total number of points in each patch

%temp_1 = zeros(N,N);
temp_2 = zeros(N,M);
temp_3 = zeros(M,N);

temp_4 = zeros(M,M);

for i=1:M
 temp_3(i,patches{i,1}')=1;   
 temp_4(i,i)=length(patches{i,1});
end

temp_1=sparse(diag(sum(temp_3)));
temp_3=-1*temp_3;
temp_2=temp_3';

L=[temp_1 temp_2;temp_3 temp_4];
L=sparse(L);




%% D matrix





%% B matrix 

for i = 1 : M
    
    labels = patches{i,1};
    tempD = zeros(dim,dim);
    tempB = zeros(dim,N+M);
    for l = 1 : length(labels)
        k = labels(l);
        e = zeros(N+M,1); e(k,1) = 1; e(N+i,1) = -1;
        x = patches{i,2}(:,l);
        tempB = tempB + (x*e');
        tempD = tempD + (x*x');
    end
    D((i-1)*dim+1 : i*dim, (i-1)*dim+1 : i*dim) = tempD;
    B((i-1)*dim+1 : i*dim,                   :) = tempB;
end




end