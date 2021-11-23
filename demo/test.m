clc;clear;close all;
%% init
addpath('./flann/');
addpath('./data/');
addpath('./ulity/');
%% prepare data
load p(Bun0);
shape= data';
scan=shape(:,1);
res= 10;
scannum=length(shape);
s=1000;
jj=1;
R_init=GrtR{1,2};
for i=1:scannum
    scan{i,1}=s*scan{i,1};              % Scale the data
end
for i=1:scannum
    %每10个取一个点
    model{i}=scan{i}(1:3,1:res:end);
    model{i} = awgn(model{i},25*jj);
end

opt.downSample=0;
opt.radii = 2000*0.002;
opt.H=10;
opt.k=5;
opt.R=eye(3);
opt.Xi=2;%控制局部向量的权重
opt.beta=1;%控制权重矩阵
opt.alphamax=30;%控制点到面权重
opt.t=[1,1,1]';

tarData=model{1,1}';
srcData=model{1,2}';
% tarData=tarData(1:4010,:);
% [T]=local_gmm(tarData,srcData,opt);

[tar_n,tar_curvature,tar_localVec,tar_localDist]=findPointNormals(tarData,opt.k);
TData=transform(srcData,opt.R,opt.t);
[T_n,T_curvature,T_localVec,T_localDist]=findPointNormals(TData,opt.k);
% [R,t]=pointToPlaneMetric(TData,tarData,tar_n);
W=ones(4026,4010);
[R,t]=pointToPlaneW(TData,tarData,tar_n,W);
result=R_init*R;