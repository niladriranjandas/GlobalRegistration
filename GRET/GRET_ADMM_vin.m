function  [X,trans,O]  = GRET_ADMM_vin( patches, N,L,D,B)


disp('ADMM running ....');

M   = size(patches,1);
dim = size(patches{1,2},1);

if(nargin<2)
[L,D,B]=createLDBMatrix2(patches,N);
end

C  = D -  B * inverse(L) * B';
disp('Created C Matrix');


% ADMM
Z=eye(M*dim);          %intialize Z with identity
Y=zeros(M*dim,M*dim); %Mdual is the dual variable
rho=.409    ;  % decreasing ro gives faster convergence but bad registration
d=3;
numIter= 250;
obj=[];
% ------------------------------------------
% The variable are 
% G - PSD
% Z - Block diagonal
% M - dual variable
% 
% rho - penality parameter
% numIter - number of iterations to be run
%-------------------------------------------

% -----------------------
% GRET - SDP using ADMM
%------------------------


for Iter=1:numIter


 G0  =  Z - (1/rho) * (C - Y);
 
   [Q,Sigma]=eig(G0);
    Sigma(Sigma<=0) = 0;
    G  = Q*Sigma*Q';
 obj=[obj,trace(C*G)];
 Z  = G - (1/rho) * Y;
 
    for i = 1:M
        Z(d*(i-1)+1:i*d, d*(i-1)+1:i*d) = eye(d);
    end
 
 
 
 
 Y= Y+rho*(Z-G);

 

end





G=(G+G')/2;

figure, axis square, semilogy(obj);
title('Objective value using ADMM vs iteration');
xlabel('iteration')
ylabel('objective value')
%determininstic rounding
[V, S] = eig(G); S=real(S); V=real(V);
Q = V(:, M*dim - dim + 1 : M*dim) * ...
    sqrt( S( M*dim - dim + 1 : M*dim, M*dim - dim + 1 : M*dim)); % top eigenvectors

% (polar decomposition) rounding to O(dim)
reconOrtho = zeros(dim,M*dim);
for i = 1 : M
    [U, ~, V]  = svd( Q(dim*(i - 1) + 1 : i*dim, :) );
    O{i}=V*U';
    reconOrtho(:, dim*(i - 1) + 1 : i*dim) = V * U';
end

% Extract coordinates
X = reconOrtho * B * inverse(L);
trans=X(:,N+1:end);

X = X(:, 1:N);


end