clc;clear;close all;

addpath('./flann/');
addpath('./data/');
addpath('./ulity/');
load p(Bun0);
shape= data';
scan=shape(:,1);
res= 10;
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
opt.downSample=0;
opt.radii = 0.02;
opt.H=10;
opt.k=5;
opt.s=1000;
[vecs,lambda,n]=computeDescriper(model{1},opt);
%% 验证下vecs
KDtree= createns(model{1}');
[corr,TD] = knnsearch(KDtree,model{1}',"k",6);
n=[1,3,4,4];
dist=norm(vecs{1,n(1)}(:,n(2)));
dist2=norm(vecs{1,n(3)}(:,n(4)));
distTD=TD(n(1),n(2)+1);
distTD2=TD(n(3),n(4)+1);
diff1=dist-distTD;
diff2=dist2-distTD2;
%% 
