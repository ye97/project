%==========================================================================
% Solve the following minimization problem:
%       min_{R, T} sum(|dot(R*p+T-q,nv)|^2)
%
% p, q, nv are all N-by-3 matrix, and nv is the unit normal at q
%
% Here the problem is solved by linear approximation to the rotation matrix
% when the angle is small.
%==========================================================================

function [R, T] = pointToPlaneMetric(p, q, nv)
% Set up the linear system
cn = [cross(p,nv,2),nv];
C = cn'*cn;
qp = q-p;
b =  [sum(sum(qp.*repmat(cn(:,1),1,3).*nv, 2));
    sum(sum(qp.*repmat(cn(:,2),1,3).*nv, 2));
    sum(sum(qp.*repmat(cn(:,3),1,3).*nv, 2));
    sum(sum(qp.*repmat(cn(:,4),1,3).*nv, 2));
    sum(sum(qp.*repmat(cn(:,5),1,3).*nv, 2));
    sum(sum(qp.*repmat(cn(:,6),1,3).*nv, 2))];
% X is [alpha, beta, gamma, Tx, Ty, Tz]
X = C\b;

cx = cos(X(1));
cy = cos(X(2));
cz = cos(X(3));
sx = sin(X(1));
sy = sin(X(2));
sz = sin(X(3));

R = [cy*cz, sx*sy*cz-cx*sz, cx*sy*cz+sx*sz;
    cy*sz, cx*cz+sx*sy*sz, cx*sy*sz-sx*cz;
    -sy,          sx*cy,          cx*cy];

T = X(4:6);
end