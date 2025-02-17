function  [Xcord,rot,trans]  = GRET_SDP_vin( cloud, N,L,D,B)


disp('SDP running....');

M   = size(cloud,1);
dim = size(cloud{1,2},1);

if(nargin<2)
[L,D,B]=createLDBMatrix2(cloud,N);
end

C  = D -  B *inverse(L)* B';
disp('Created C Matrix');

%call SDP-CVX
G = zeros(M*dim,M*dim);
tic;
cvx_precision best;
I = eye(dim);
cvx_begin sdp                             
variable G(M*dim, M*dim) symmetric
minimize( trace(C*G))
subject to
G >= 0 ;
for i = 1 : M
    G(dim*(i - 1) + 1 : dim*i, dim*(i - 1) + 1 : dim*i) == I;  
end;
cvx_end
fprintf('Finished computing the SDP');

disp('^^^^^^^^^^^^^^^^^');
disp(trace(C*G));
disp('^^^^^^^^^^^^^^^^^');

% determininstic rounding
[V, S] = eig(G); S=real(S); V=real(V);


Q = V(:, M*dim - dim + 1 : M*dim) * ...
    sqrt( S( M*dim - dim + 1 : M*dim, M*dim - dim + 1 : M*dim)); % top eigenvectors

reconOrtho = zeros(dim,M*dim);
for i = 1 : M
    [U, ~, V]  = svd( Q(dim*(i - 1) + 1 : i*dim, :) );
    reconOrtho(:, dim*(i - 1) + 1 : i*dim) = V * U';
    rot{i}=V*U';
end

% Extract coordinates
Xcord = reconOrtho * B * inverse(L);
Xcord = Xcord(:, 1:N);

%Extract Translations
trans=Xcord(:,N+1:end);
end

