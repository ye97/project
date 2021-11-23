function [R, T] = pointToPlaneMetric(p, q, nv,W)
%==========================================================================
% Solve the following minimization problem:
%       min_{R, T} sum(W|dot(R*p+T-q,nv)|^2)
%
% p, q, nv are all N-by-3 matrix, and nv is the unit normal at q
%
% Here the problem is solved by linear approximation to the rotation matrix
% when the angle is small.
%==========================================================================
[N,D]=size(q);
[M,D]=size(p);
% N=4026,M=4010;

% cn = [cross(nv,p,2),nv];
% cn=bsxfun(@cross,nv,p);
en=q.*nv;
nv=repmat(nv,1,M);
nv=reshape(nv,N,3,M);
nv=permute(nv,[1,3,2]);
p=repmat(p,1,N);
p=reshape(p, M,3,N);
p=permute(p,[3,1,2]);
% cn=cross(p,nv);
cn=bsxfun(@cross,p,nv);
% for i=1:N
%     for j=1:M
%         TEMP(i,j,:)=cross(p(i,j,:),nv(i,j,:));
%     end
% end
% result=TEMP-cn;
% sum(sum(result));
% disp(p);
