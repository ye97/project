N=4;
M=5;
TData=rand(M,3);
ORIGIN_T=TData;
TData=repmat(TData,1,N);
TData=reshape(TData,[M,3,N]);
TData=permute(TData,[3,1,2]);
a=[1,2,3];
a=a.^(1/2);
p=rand(3,3);
q=rand(4,3);
temp=p-q;

A=rand(N,M,3);
B=rand(N,M,3);
TEST= bsxfun(@)