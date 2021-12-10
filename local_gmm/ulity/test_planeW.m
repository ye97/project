function [R, T,sigma2] = test_planeW(p, q, nv,W)
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


origin_W=W;
W=permute(W,[2,1]);
W=reshape(W,N*M,1);

W=W.^(1/2);

en=sum(q.*nv,2);
en=repmat(en,1,M);
en=permute(en,[2,1]);
en=reshape(en,N*M,1);
en_NM=en;
B=-en.*W;

d=repelem(nv,M,1).*repmat(p,N,1);
d=sum(d,2);
d=d.*W;
B=B+d;

W=repmat(W,1,3);
A=cross(repelem(nv,M,1),repmat(p,N,1));
A=A.*W;

origin_nv=nv;

nv=repmat(nv,1,M);
nv=reshape(nv,[N,3,M]);
nv=permute(nv,[2,3,1]);
nv=reshape(nv,[3,N*M]);
nv=permute(nv,[2,1]);

nv=-nv.*W;


A=[A  nv];
X=A\B;
% for i=1:N
%     for j=1:M
%         TEMP(i,j,:)=cross(p(i,j,:),nv(i,j,:));
%     end
% end
% result=TEMP-cn;
% sum(sum(result));

cx = cos(X(1));
cy = cos(X(2));
cz = cos(X(3));
sx = sin(X(1));
sy = sin(X(2));
sz = sin(X(3));
R=[1,0,0;0,cx,sx;0,-sx,cx;]*[cy,0,-sy;0,1,0; sy,0,cy;]*[cz,sz,0; -sz,cz,0;0,0,1;];
% R = [cy*cz, sx*sy*cz-cx*sz, cx*sy*cz+sx*sz;
%     cy*sz, cx*cz+sx*sy*sz, cx*sy*sz-sx*cz;
%     -sy,          sx*cy,          cx*cy];


T = X(4:6)';
sigma2=sum((A*X+en_NM).^2)/sum(sum(origin_W));



