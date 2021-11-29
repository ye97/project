N=3;
M=4;
% z=rand(N,M,6);
% a=permute(z,[3,2,1]);
% % result=reshape(a,N*M,6);
% result=reshape(a,6,N*M);
% result=result';
% 
% 
% % N*3->N*M*3
% m=rand(N,3);
% origin_m=m;
% m=repmat(m,1,M);
% m=reshape(m,N,3,M);
% m=permute(m,[1,3,2]);
% 
% %M*3->N*M*3
% p=rand(M,3);
% p_NM=repmat(p,1,N);
% p_NM=reshape(p_NM, M,3,N);
% p_NM=permute(p_NM,[3,1,2]);


%  N*3->(N*M)*3
% 每行重复m次
nv=rand(N,3);
origin_E=nv;
nv=repmat(nv,1,M);
nv=reshape(nv,[N,3,M]);
nv=permute(nv,[2,3,1]);
nv=reshape(nv,[3,N*M]);
nv=permute(nv,[2,1]);

W=rand(N,M);
origin_w=W;
W=permute(W,[2,1]);
W=reshape(W,N*M,1);