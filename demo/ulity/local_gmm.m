function [T]=local_gmm(tarData,srcData,opt)
%% Check the input options and set the defaults

if ~isfield(opt,'H') || isempty(opt.H), opt.max_it = 10; end
if ~isfield(opt,'k') || isempty(opt.k), opt.k = 5; end
if ~isfield(opt,'downSample') || isempty(opt.downSample), opt.downSample=0; end
if ~isfield(opt,'s') || isempty(opt.s), opt.s = 1000; end
if ~isfield(opt,'radii') || isempty(opt.radii), opt.radii = 0.02; end
if ~isfield(opt,'radii') || isempty(opt.radii), opt.radii = 0.02; end
w=0.2;
max_it=200;
[tar_n,tar_curvature,tar_localVec,tar_localDist]=findPointNormals(tarData,opt.k);
[N,D]=size(tarData);
[M,D]=size(srcData);
if ~exist('sigma2','var') || isempty(sigma2) || (sigma2==0), 
    sigma2=(M*trace(tarData*tarData')+N*trace(srcData*srcData')-2*sum(srcData)*sum(tarData)')/(M*N*D);
end

iter=0;
WeightMatrix=[];
ntol=10;

%% 
while (iter<max_it) && (sigma2 > 1e-8)% && (ntol<1e-8) 
    TData=transform(srcData,opt.R,opt.t);
    [T_n,T_curvature,T_localVec,T_localDist]=findPointNormals(TData,opt.k);
    %% E step
    paiMatrix=computePai(tar_localVec,T_localVec, opt.beta);
    alpha=compute_alpha(tar_curvature,opt.alphamax);
    gloDist=compute_gloDist(tarData,TData,tar_n,tar_curvature,alpha);
    P_prior =comput_Prior(paiMatrix,gloDist,w,sigma2);
    E_con=compute_E(P_prior,alpha,gloDist,sigma2);
    [R,t]=update(TData,tarData,tar_n,P_prior,alpha,gloDist,sigma2);
    %% M step
    
    iter=iter+1;
end
T=WeightMatrix;
disp(iter);
disp(sigma2);