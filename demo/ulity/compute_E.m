function    [E_con]=compute_E(P_prior,alpha,gloDist,sigma2)
[N,M]=size(P_prior);
alpha=repmat(alpha,1,M);
E_con=sum(sum(P_prior.*alpha.*gloDist./2/sigma2))+N*M*log(sigma2)/2;