function [vecs,n,lambda,srcCloud] = computeDescriper(srcCloud,opt)
%This code is the Matlab implimentation of the paper, 
%"Fast Descriptors and Correspondence Propagation for Robust Global Point Cloud Registration,"
%IEEE transactions on Image Processing, 2017.
%This code should be used only for academic research.
%any other useage of this code should not be allowed without Author agreement.
% If you have any problem or improvement idea about this code, please
% contact with Guang JIANG, Xidian University. gjiang@mail.xidian.edu.cn.



%% parameter configuration for flann search
params.algorithm = 'kdtree';
params.trees = 8;
params.checks = 64;

srcData = srcCloud;

% 半径0.5间隔
radii = opt.radii;

if opt.downSample==1
    srcCloudDown = pcdownsample(srcCloud, 'gridAverage', gridStep);
   
    srcSeed = srcCloudDown;
   
else
    srcSeed = srcData;
    
end



%% compute descriptors for seed points in the source point cloud


% IDX = RANGESEARCH(X,Y,RADIUS) finds all the points in X that are
% within distance RADIUS for points in Y.
% 找到data中半径为r的所有点
srcIdx = rangesearch(srcData',srcSeed',radii);
% cellfun 对每个cell进行函数操作
% idxsize表示srcseed有多少个邻近点
idxSz = cellfun(@length,srcIdx,'uni',true);
srcIdx = srcIdx(idxSz>10);
% srcIdxK=srcIdx{:}(1:k);
% srcSeed至少有10邻近点
srcSeed = srcSeed(:,idxSz>opt.H);
% 种子点个数
M = sum(idxSz>opt.H);
N=sum(idxSz<=opt.H);
% 1-M的下标
idx = num2cell((1:M)');
% 匿名函数设置x，y


[lambda,n] = cellfun(@(x,y)svdCov(x,y,srcData,srcSeed),srcIdx,idx,'uni',false);
vecs=[];

for i=1:M
    vec=[0,0,0]';
    for j=2:opt.k+1
        vec_ij=srcSeed(:,i)-srcData(:,srcIdx{i,1}(1,j));
%         局部距离约束
%         vec_ij=norm(vec_ij);
%         vec=[vec,vec_ij];
          vec=vec+vec_ij./norm(vec_ij);
    end
    vecs{i}=vec;
end
vecs=vecs';
srcCloud=srcSeed;

vecs=cell2mat(vecs);
vecs=reshape(vecs,M,3,1);


lambda=cell2mat(lambda);
lambda=reshape(vecs,M,3,1);

n=cell2mat(n);
n=reshape(vecs,M,3,1);
