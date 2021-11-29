clc;clear;close all;
%% init
addpath('./flann/');
addpath('./data/');
addpath('./ulity/');
%% prepare data
load p(Bun0);
shape= data';
scan=shape(:,1);
res= 20;
scannum=length(shape);
s=1000;
jj=1;
for i=1:scannum
    scan{i,1}=s*scan{i,1};              % Scale the data
end
for i=1:scannum
    %每10个取一个点
    model{i}=scan{i}(1:3,1:res:end);
    model{i} = awgn(model{i},25*jj);
end
%% config parameters
opt.downSample=0;
opt.radii = 2000*0.002;
opt.H=10;
opt.k=5;
opt.R=eye(3);
opt.Xi=2;%控制局部向量的权重
opt.beta=1;%控制权重矩阵
opt.alphamax=1;%控制点到面权重
opt.t=[0,0,0]';

tarData=model{1,1}';
srcData=model{1,2}';
theta = pi/24;
rot   = [cos(theta) sin(theta) 0; ...
          -sin(theta) cos(theta) 0; ...
                   0          0  1];
trans = [5 5 10];
tform1 = rigid3d(rot, trans);
tarData=srcData*rot+trans;
tic;
t1=toc;
ptCloud1=pointCloud(tarData);
pcwrite(ptCloud1,'tar.ply');

ptCloud2=pointCloud(srcData);
pcwrite(ptCloud2,'src.ply');
ptCloud_src = pcread('src.ply');
ptCloud_tar = pcread('tar.ply')
[tform2, movingReg, rmse] = pcregistericp(ptCloud_src, ptCloud_tar, 'Extrapolate', true,'Metric',"pointToPlane");

[tar_n,tar_curvature,tar_localVec,tar_localDist]=findPointNormals(tarData,opt.k);
[R1,t1]=pointToPlaneMetric(srcData,tarData,tar_n);
paiMatrix=ones(2005,2005)./2005;
[R_w,t_w,sigma_w]=pointToPlaneW(srcData,tarData,tar_n,paiMatrix);
[T_cpd,~,~]=pcregistercpd(ptCloud_src,ptCloud_tar,'Transform',"rigid");

[T]=local_gmm(srcData,tarData,opt);
% 
% W2=ones(4026,4026);

t2=toc;
time=t2-t1;

% [tar_vecs,tar_n,tar_lambda,tarData]=computeDescriper(tarData,opt);
% 
% shape=transform(srcData,opt.R,opt.t);
% 
% [src_vecs,src_n,src_lambda,shape]=computeDescriper(shape,opt);

% WeightMatrix=computeWeightMatrix(tar_vecs,src_vecs,tar_n,src_n,opt.beta);
% % 查看最大权重组件
% maxN=max(WeightMatrix);
% maxNM=max(maxN);


%% 验证下vecs
% KDtree= createns(model{1}');
% [corr,TD] = knnsearch(KDtree,model{1}',"k",6);
% n=[1,3,4,4];
% dist=norm(vecs{1,n(1)}(:,n(2)));
% dist2=norm(vecs{1,n(3)}(:,n(4)));
% distTD=TD(n(1),n(2)+1);
% distTD2=TD(n(3),n(4)+1);
% diff1=dist-distTD;
% diff2=dist2-distTD2;
%% 
