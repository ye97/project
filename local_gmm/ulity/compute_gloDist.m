function [gloDist]=compute_gloDist(tarData,TData,tar_n,alpha)

[N,D]=size(tarData);
[M,D]=size(TData);
tar_n=repmat(tar_n,1,M);
tar_n=reshape(tar_n,[N,3,M]);
tar_n=permute(tar_n,[1,3,2]);
% for i=1:N
%     for j=1:M
%         DIST(i,j)=norm(tarData(i,:)-srcData(j,:));
%     end
% end
TData=repmat(TData,1,N);
TData=reshape(TData,[M,3,N]);
TData=permute(TData,[3,1,2]);

tarData=repmat(tarData,1,M);
tarData=reshape(tarData,[N,3,M]);
tarData=permute(tarData,[1,3,2]);

gloDist=(tarData-TData).*tar_n;
gloDist=sum(gloDist,3);
% gloDist=gloDist.^2;
alpha=repmat(alpha,1,M);
gloDist=gloDist.*alpha;
% gloDist=sqrt(sum((tarData-srcData).^2,3));
% gloDist=squeeze(gloDist);
% alpha=repmat(alpha,1,M);
% gloDist=gloDist.^2.*alpha;


% 传统距离
% gloDist=sqrt(sum((tarData-srcData).^2,3));
% test=DIST-gloDist;

