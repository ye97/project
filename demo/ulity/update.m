function [R,t]=update(TData,tarData,tar_n,P_prior,alpha,gloDist,sigma2)
[N,M]=size(P_prior);
alpha=repmat(alpha,1,M);
P_prior=power(P_prior,1/2);
alpha=power(alpha,1/2);
TData=TData.*P_prior.*alpha;
[R,t]=pointToPlaneMetric(TData, tarData, tar_n);

