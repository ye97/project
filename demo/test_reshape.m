N=3;
M=4;
z=rand(N,M,6);
a=permute(z,[3,2,1]);
result=reshape(a,6,N*M);
result=result';