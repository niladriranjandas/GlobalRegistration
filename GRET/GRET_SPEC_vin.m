function  [X,rot,trans]  = GRET_SPEC_vin( patches, N,L,D,B)


disp('SPEC running ....');

M   = size(patches,1);
dim = size(patches{1,2},1);

if(nargin<2)
[L,D,B]=createLDBMatrix2(patches,N);
end

C  = D -  B * inverse(L) * B';

disp('Created C Matrix');

[V, S] = eigs(C,90); S=real(S); V=real(V);
  
  V=V(:,M*dim-dim+1:end);
  W=sqrt(M)*(V)';

for i=1:M
    [ui,~,vi]=svd(W(:,(i-1)*dim+1 : i*dim));
    O{i}=ui*vi';
    rot{i}=O{i};

end


O=horzcat(O{:});

% Extract coordinates
X = O * B * inverse(L);
trans=X(:,N+1:end); %Extract Translation
disp(trans)

X = X(:, 1:N);



end
