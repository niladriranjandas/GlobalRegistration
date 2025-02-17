function [J,B,D] = const_mat(Patches,N)
    %% Function Description
    % This is function is used to create the constant matrices J,B,D for
    % solving GRET-SDP or GRET-SPEC
    % Inputs:-
    % Patches => Cell: 2; Cell-1: Identity of node; Cell-2: Local   
    %            co-ordinates
    % N       => Number of sensors
    % Output:-
    % J: N+M-by-N+M matrix
    % B:  Md-by-N+M matrix
    % D:  md-by-Md  matrix
    %% Parameters
    M = size(Patches,1);                     %Calculating number of patches
    d = size(Patches{1,2},1);                                    %Dimension
    %% Problem initialization
    J = zeros(N+M);
    B = zeros(M*d,(N+M));
    D = zeros(M*d);
    %% Updating the variables
    for i = 1:M
        node = Patches{i,1};                %Nodes that is in the ith patch 
        for j = 1:length(node)
            % Resetting the variables
            eij = zeros(N+M);
            B1 = zeros(M*d,(N+M));
            D1 = zeros(M*d);
            % C Formation
            k = node(j);                                      % Node Number
            xki = Patches{i,2}(:,j);                      % Position Vector
            xkii = xki*xki';
            % Updating eij, B1, D1- Sensors
            eij(k,k) = 1;
            eij(N+i,N+i) = 1;
            eij(k,N+i) = -1;
            eij(N+i,k) = -1;
            B1(((i-1)*d+1):(i*d),k) = xki;
            B1(((i-1)*d+1):(i*d),(N+i)) = -xki;
            D1(((i-1)*d+1):(i*d),((i-1)*d+1):(i*d)) = xkii;
            % Main updates for Sensor Node- J, B, D
            J = J+eij;
            B = B+B1;
            D = D+D1;
        end
    end
    J = (J+J')/2;
    D = (D+D')/2;
end