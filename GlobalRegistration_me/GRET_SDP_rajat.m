function [Z,O,phi,Err] = GRET_SDP(J,B,D,X)
    % Function describtion
    % Inputs:-
    % Patches => Cell: 2; Cell-1: Identity of node; Cell-2: Local   
    %            co-ordinates
    % N       => Number of sensors
    % X       => d-by-N matrix where ith column represents the original 
    %            position of the ith sensor 
    % Output:-
    % Z: d-by-(N+M) matrix where first N columns represent the the
    %    recovered positions of the sensors and (N+i)th column represents
    %    the ith patch translation vector
    % O: d-by-Md matrix where O=[O1 ... OM]: Oi represents the ith patch
    %    rotation matrix
    % phi: Objective function value
    % Err: Registration error
    %% Parameters
    M = size(J,1);
    [d,N] = size(X);                        % d: Dimension; N: No. of nodes
    M = M-N;                                % M: No. of patches
    %% Creating C matrix
    C = D-B*pinv(J)*B';
    C  = (C+C')/2;                             % Making it symmetric if not
    %% Calling CVX Package
    cvx_precision best;
    cvx_begin
    variable G(M*d,M*d) semidefinite                   % Defining variables
    minimize(trace(C*G))                               % Objective function
    subject to
    % Constraints
    for i = 0:(M-1)
        G(i*d+(1:d),i*d+(1:d)) == eye(d);
    end
    cvx_end
    %% O(d) Recovery
    % Note: Due to numerical precision, if G is not comming symmetric,
    % because of that we are doing-
    G = (G+G')/2;
    rk = rank(G,0.0001);
    fprintf('Rank of G matrix is: %d \n',rk)
    [V,e] = eigs(G,d,'LA');
    W = sqrt(e)*V';
    O = zeros(d,M*d);
    for i = 0:(M-1)
        [U2,~,V2] = svd(W(1:d,i*d+(1:d)));
        O(1:d,i*d+(1:d)) = U2*V2';
    end
    %% Final Result Computation
    Z = O*B*pinv(J);
    phi = trace(Z*J*Z'-2*O*B*Z'+O*D*O');
    X_est = Z(:,1:N);
    [~,Err,~,~] = Globalregistration_Arun(X_est,X,0);
end