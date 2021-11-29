function [Prior]=comput_Prior(paiMatrix,gloDist,w,sigma2)
[N,M]=size(gloDist);
Prior=zeros(N,M);
top=paiMatrix.*exp(-gloDist/(2*sigma2));
c  = power(2*pi*sigma2, 1/2) *( w /(1-w) /N );

pMatrix_col=sum(top,2)+c;
pMatrix_col=repmat(pMatrix_col,1,M);
% pMatrix_bot=permute(pMatrix_col,[2,1]);
Prior=top./(pMatrix_col);
