function [gloDist]=compute_gloDist(tarData,srcData,tar_n,tar_curvature,alpha)
[N,D]=size(tarData);
[M,D]=size(srcData);
% for i=1:N
%     for j=1:M
%         DIST(i,j)=norm(tarData(i,:)-srcData(j,:));
%     end
% end
srcData=repmat(srcData,1,1,N);
srcData=permute(srcData,[3,1,2]);
tarData=repmat(tarData,1,1,M);
tarData=permute(tarData,[1,3,2]);
gloDist=sum((tarData-srcData).^2,3);
% gloDist=sqrt(sum((tarData-srcData).^2,3));
gloDist=squeeze(gloDist);
alpha=repmat(alpha,1,M);
gloDist=gloDist.*alpha;
% 传统距离
% gloDist=sqrt(sum((tarData-srcData).^2,3));
% test=DIST-gloDist;

