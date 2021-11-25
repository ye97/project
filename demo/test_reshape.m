N=3;
M=4;
z=rand(N,M,6);
a=permute(z,[3,2,1]);
result=reshape(a,6,N*M);
result=result';


% N*3->N*M*3
m=rand(N,3);
origin_m=m;
m=repmat(m,1,M);
m=reshape(m,N,3,M);
m=permute(m,[1,3,2]);