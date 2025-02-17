function [R, T] = ...
    MVSE3EstimateL2Krishnan(MInit, pts, ijMap, corrMap, maxIters, tol)
%MVSE3EstimateL2Krishnan estimate the global motion matricies
%    Implements [1] to estimate the global motion matricies (MglobalEst) from
%    a set of given point clouds (pts), the view graph edge set (ijMap),
%    correspondences (corrMap).
%    
%    [MglobalEst, viewsAligned, err, errs, A, B, C] =
%    MVSE3ESTIMATEL2KRISHNAN(MInit, pts, ijMap, corrMap, maxIters, tol)
%    Inputs:
%    MInit: Matrix of size (4 x 4 x nViews) Initial global motion estimate.
%    pts: cell of size (1 x nViews) with each cell containing a matrix of
%         size (3 x n_i) where n_i is the number of points in the point cloud
%    ijMap: matrix of size (nRels x 2) where nRels is the number of relations
%    corrMap: cell of size (1 x nRels). If [i,j] = ijMap(k,:), corrMap(k,1)
%    points of pts{i} corresponds to corrMap(k,2) points of pts{j}.
%    maxIters: maximum number of newton steps.
%    convergence criterion.
%    Outputs:
%    MglobalEst: (4 x 4 x nViews) global motion estimate
%    viewsAligned: aligned pts
%    err: final error
%    errs: error at each iteration
%    A, B, C: see [1] for definitions.
%
%    [1] Shankar Krishnan, Pei Yean Lee, John B Moore, and Suresh
%    Venkatasubramanian. Global registration of multiple 3d point sets via
%    optimization-on-a-manifold. In Symposium on Geometry Processing, pages
%    187–196, 2005
%
%    Copyright (C) 2014 Anil C R <cr.anil@gmail.com>
%
%    This library is free software; you can redistribute it and/or
%    modify it under the terms of the GNU General Public License as published
%    by the Free Software Foundation; either version 2.1 of the License, or
%    (at your option) any later version.
%
%    This library is distributed in the hope that it will be useful, but
%    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
%    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
%    for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program. If not, see <http://www.gnu.org/licenses/>.
%
%    Authored by Anil C R <cr.anil@gmail.com>

nViews = length(pts);
alpha = 0.01;
beta = 0.5;
nRels = size(ijMap,1);

A = zeros(3*nViews, 3*nViews);
B = zeros(3*nViews,nViews);
C = zeros(nViews, nViews);

Qx = [ 0  0  0
       0  0 -1
       0  1  0];
Qy = [ 0  0  1
       0  0  0
       -1  0  0];
Qz = [ 0 -1  0
       1  0  0
       0  0  0];
    function R = update(wopt, R)
    for iteration=1:nViews
        w = wopt(3*iteration-2:3*iteration);
        R(:,3*iteration-2:3*iteration) = R(:,3*iteration-2:3*iteration)*expm(w(1)*Qx+w(2)*Qy+w(3)*Qz);
    end
    end

N = 0;

for edg=1:nRels
    corr = corrMap{edg};
    nCorrs = size(corr,1);
    i = ijMap(edg,1);
    j = ijMap(edg,2);
    vij = pts{i}(:, corr(:,1));
    vji = pts{j}(:, corr(:,2));
    A(3*i-2:3*i,3*i-2:3*i) = A(3*i-2:3*i,3*i-2:3*i) + vij*vij';
    A(3*i-2:3*i,3*j-2:3*j) = A(3*i-2:3*i,3*j-2:3*j) - vij*vji';
    A(3*j-2:3*j,3*i-2:3*i) = A(3*j-2:3*j,3*i-2:3*i) - vji*vij';
    A(3*j-2:3*j,3*j-2:3*j) = A(3*j-2:3*j,3*j-2:3*j) + vji*vji';
    B(3*i-2:3*i,i) = B(3*i-2:3*i,i) + sum(vij,2);
    B(3*j-2:3*j,j) = B(3*j-2:3*j,j) + sum(vji,2);
    B(3*j-2:3*j,i) = B(3*j-2:3*j,i) - sum(vji,2);
    B(3*i-2:3*i,j) = B(3*i-2:3*i,j) - sum(vij,2);
    C(i,i) = C(i,i) + nCorrs;
    C(j,j) = C(j,j) + nCorrs;
    C(i,j) = -nCorrs;
    C(j,i) = -nCorrs;
    N = N + nCorrs;
end

RInit = zeros(3,3*nViews);
TInit = zeros(3,nViews);

for i=1:nViews
    RInit(:,3*i-2:3*i) = MInit(1:3,1:3,i);
    TInit(:,i) = MInit(1:3,4,i);
end

fp = trace(RInit*A*RInit'+2*RInit*B*TInit'+TInit*C*TInit')/N;
errs = zeros(1, maxIters+1);
errs(1) = fp;
fprintf('\t\tIteration: %d, ERROR: %05.016f\n', 0, fp);

B = B(:,2:end);
C = C(2:end,2:end);
M = A-B*(C\B');
Qe = cell(nViews,1);
for i=1:nViews
    ei = zeros(nViews,1);
    ei(i) = 1;
    Qe{i} = [kron(ei,Qx);kron(ei,Qy);kron(ei,Qz)];
end
Q = blkdiag(Qe{:});
R = RInit;
for iter=1:maxIters
    J = kron(R,eye(3*nViews))*Q;
    g = 2*J'*vec(M*R');
    H = J'*kron(eye(3),M)*J;
    %   Hr = -Q'*kron(eye(3*nViews-3), M*R'*R)*Q;
    wopt = -H\g;
    
    %   backtrack
    t = 1;
    fp = trace(R*M*R');
    for backiter=1:100
        Rp = update(t*wopt,R);
        fnow = trace(Rp*M*Rp');
        if fnow <= fp + t*alpha*g'*wopt
            break;
        else
            t = beta*t;
        end
    end
    wopt = t*(wopt);
    R = update(wopt,R);
    T = -R*(B/C);
    size(R)
    size(M)
    
    
    
    size(B)
    size(T)
    size(C)
    fp = trace(R*M*R'+2*R*B*T'+T*C*T')/N;
    T=[zeros(3,1) T];
    errs(iter+1) = fp;
    
    cnd = norm(wopt);
    fprintf(['\t\tIteration: %d, ERROR: %05.016f,',...
             ' norm(W_inc): %05.016f, norm(g): %05.016f\n'],...
            iter, fp, norm(wopt), norm(g));
    if cnd < tol, break;end
end
% errs = errs(1:iter+1);
% MglobalEst = zeros(4,4,nViews);
% MglobalEst(:,:,1) = eye(4);
% for i=1:nViews
%     MglobalEst(1:3,1:3,i) = R(1:3,1:3)\R(:,3*i-2:3*i);
%     MglobalEst(1:3,4,i) = R(1:3,1:3)\T(3*i-2:3*i)';
%     MglobalEst(4,4,i)=1;
% end

% viewsAligned = cell(size(pts));

% for k = 1:size(MglobalEst,3)
%     viewsAligned{k} = transform(pts{k}, MglobalEst(:,:,k));
% end

% err = calcErrorl2(MglobalEst,pts,ijMap,corrMap);
end