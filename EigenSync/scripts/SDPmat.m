
function X= SDPmat(X0)


[m,n]=size(X0);


cvx_begin sdp quiet
    variable X(n,n) symmetric
    X == semidefinite(n);
    maximize(trace(X*X0))
    diag(X)==1 
    
cvx_end



end