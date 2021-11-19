clc;clear;close all;
addpath('./flann/');
addpath('./data/');
addpath('./ulity/');
load p(Bun0);
scannum=10;
scan=data';
res=10;
for i=1:scannum
    scan{i,1}=s*scan{i,1};              % Scale the data
end
for i=1:scannum
    %每10个取一个点
    model{i}=scan{i}(1:3,1:res:end);
end
tarData=model{1,1};
srcData=model{1,2};
opt.downSample=0;
opt.radii = 2000*0.002;
opt.H=10;
opt.k=5;
opt.s=1000;
opt.R=eye(3);
opt.beta=1;%控制权重矩阵
opt.alpha=30;%控制点到面权重
opt.t=[1,1,1]';
tarData=model{1,1};
srcData=model{1,2};


[DescriptorTar,tarData]=computeDescriper(tarData,opt);
% martrixN=cell2mat(DescriptorTar{:,2});
% ny{:,:}=DescriptorTar{:,2}{:,:};
d=cell2mat(DescriptorTar);

% [R,t]=pointToPlaneMetric(srcData',tarData',DescriptorTar{:,2});