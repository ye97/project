function [R, T,sigma2] = pointToPlaneMetric(p, q, nv,W)
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
sigma2=0;
% cn = [cross(nv,p,2),nv];
% cn=bsxfun(@cross,nv,p);

origin_W=W;
W=repmat(W,1,1,3);
W=W.^1/2;

en=sum(q.*nv,2);
en_NM=repmat(en,1,M);

en_NM=reshape(en_NM,N,1,M);
en_NM=permute(en_NM,[1,3,2]);
en_NM=repmat(en_NM,1,1,3);


nv_NM=repmat(nv,1,M);
nv_NM=reshape(nv_NM,N,3,M);
nv_NM=permute(nv_NM,[1,3,2]);


p_NM=repmat(p,1,N);
p_NM=reshape(p_NM, M,3,N);
p_NM=permute(p_NM,[3,1,2]);
% cn=cross(p,nv);
cn=bsxfun(@cross,p_NM,nv_NM);

cn=cn.* W;%前三列
nv_NM=-nv_NM.*W;

A=cat(3,cn,nv_NM);

% for i=1:N
%     for j=1:M
%         TEMP((i-1)*M+j,:)=A(i,j,:);
%     end
% end

A=reshape(A,[M*N,6]);
% result=TEMP-cn;
sum(sum(result));
B=repmat(en,1,M);
B=-B.*origin_W;
B=reshape(B,[M*N,1]);
% 
% % W=reshape(W,[M*N,1]);

X=inv((A'*A))*A'*B;
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

R = [cy*cz, sx*sy*cz-cx*sz, cx*sy*cz+sx*sz;
    cy*sz, cx*cz+sx*sy*sz, cx*sy*sz-sx*cz;
    -sy,          sx*cy,          cx*cy];

T = X(4:6);

x=X(1:3);
x=repmat(x,1,M*N);
x=reshape(x,3,N,M);
x=permute(x,[2,3,1]);


t=repmat(T,1,M*N);
t=reshape(t,3,N,M);
t=permute(t,[2,3,1]);


sigma2=sum(sum((cn.*x+nv_NM.*t+en_NM).^2))/M*N;
sigma2=sigma2(:,:,1);
