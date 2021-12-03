function [gloDist]=compute_gloDist(tarData,TData,tar_n,alpha)

[N,D]=size(tarData);
[M,D]=size(TData);
%% 点到面距离
y=tarData*tar_n';
y_dist=diag(y);
x_dist=tar_n*TData';
y_dist=repmat(y_dist,[1,M]);
gloDist=-x_dist+y_dist;
gloDist=gloDist.^2;
% alpha=repmat(alpha,1,M);
% gloDist=gloDist.*alpha;
%% 
% tarData=repmat(tarData,M,1);
% TData=repmat(TData,N,1);
% gloDist=sqrt(sum((tarData-TData).^2,3));
% gloDist=squeeze(gloDist);
% % alpha=repmat(alpha,1,M);
% gloDist=gloDist.^2.*alpha;


% 传统距离
% gloDist=sqrt(sum((tarData-srcData).^2,3));
% test=DIST-gloDist;

