function [T]=local_gmm(tarData,srcData,opt)
%% Check the input options and set the defaults

if ~isfield(opt,'H') || isempty(opt.H), opt.max_it = 10; end
if ~isfield(opt,'k') || isempty(opt.k), opt.k = 5; end
if ~isfield(opt,'downSample') || isempty(opt.downSample), opt.downSample=0; end
if ~isfield(opt,'s') || isempty(opt.s), opt.s = 1000; end
if ~isfield(opt,'radii') || isempty(opt.radii), opt.radii = 0.02; end
if ~isfield(opt,'radii') || isempty(opt.radii), opt.radii = 0.02; end

max_it=200;
[tar_n,tar_curvature,tar_localVec,tar_Dist]=findPointNormals(tarData,opt.k);
[N,D]=size(tarData);
[M,D]=size(srcData);
if ~exist('sigma2','var') || isempty(sigma2) || (sigma2==0), 
    sigma2=(M*trace(tarData*tarData')+N*trace(srcData*srcData')-2*sum(srcData)*sum(tarData)')/(M*N*D);
end





iter=0;
WeightMatirx=[];

%%
while (iter<max_it) && (sigma2 > 1e-8) 
    TData=transform(srcData,opt.R,opt.t);
    [src_n,src_curvature,src_localVec,src_Dist]=findPointNormals(TData,opt.k);
    
    WeightMatrix=computeWeightMatrix(tar_localVec,src_localVec,opt.beta);
    P_prior = WeightMatrix;
    iter=iter+1;
end
T=WeightMatrix;