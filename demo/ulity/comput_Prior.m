function [Prior]=comput_Prior(paiMatrix,gloDist,w,sigma2)
[N,M]=size(gloDist);
Prior=zeros(N,M);

top=paiMatrix.*exp(-gloDist./(2*sigma2));
c  = power(2*pi*sigma2, 1/2) *( w /(1-w) /N );
pMatrix_col=sum(paiMatrix,1);
pMatrix_col=repmat(M,1);
pMatrix_bot=permute(pMatrix_col,[2,1]);
Prior=top./(pMatrix_bot+c);
