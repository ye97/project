clc;clear;close all;
%% init
addpath('./data/');
addpath("./ulity/");
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
%     model{i} = awgn(model{i},25*jj);
end

srcData=model{1,2}';
theta = pi/12;
rot   = [cos(theta) sin(theta) 0; ...
          -sin(theta) cos(theta) 0; ...
                   0          0  1];
trans = [10 5 10];
tform1 = rigid3d(rot, trans);
tarData=srcData*rot+trans;


%% config parameters
opt.downSample=0;
opt.radii = 2000*0.002;
opt.H=10;
opt.k=5;
opt.R=eye(3);
opt.Xi=1.2;%控制局部向量的权重
opt.beta=3;%控制权重矩阵
opt.alphamax=1;%控制点到面权重
opt.t=[0,0,0]';
opt.R=eye(3);
opt.normalize=1;

%% local cpd
tic;
time1=toc;
cpd_plot_iter(tarData,srcData);
hold off;

[tar_n,tar_curvature,tar_localVec,tar_localDist]=findPointNormals(tarData,opt.H,opt.k);
[N,D]=size(tarData);
[M,D]=size(srcData);
paiMatrix=ones(N,M)./M;
[opt.R,opt.t,sigma]=test_planeW(srcData,tarData,tar_n,paiMatrix);
T=local_cpd(srcData,tarData,opt);
time2=toc;
time=time2-time1;






