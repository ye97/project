function [Prior, L]=comput_Prior(paiMatrix,gloDist,w,sigma2)
[N,M]=size(gloDist);
D=3;
Prior=zeros(N,M);
top=paiMatrix.*exp(-gloDist/(2*sigma2));
c  = power(2*pi*sigma2, 1/2) *( w /(1-w) /N );

sum_temp=sum(top,2)+c;

% pMatrix_bot=permute(pMatrix_col,[2,1]);
Prior=top./repmat(sum_temp,1,N);
Prior=Prior';
Np=sum(sum(Prior));
% L =  sum(log(sum_temp)) + D*N*log(sigma2)/2;
L=  sum(sum(Prior.*gloDist))/2*sigma2+ Np*log(sigma2)/2;