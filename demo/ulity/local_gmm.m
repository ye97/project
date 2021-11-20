function [T]=local_gmm(tarData,srcData,opt)
%% Check the input options and set the defaults

if ~isfield(opt,'H') || isempty(opt.H), opt.max_it = 10; end
if ~isfield(opt,'k') || isempty(opt.k), opt.k = 5; end
if ~isfield(opt,'downSample') || isempty(opt.downSample), opt.downSample=0; end
if ~isfield(opt,'s') || isempty(opt.s), opt.s = 1000; end
if ~isfield(opt,'radii') || isempty(opt.radii), opt.radii = 0.02; end
if ~isfield(opt,'radii') || isempty(opt.radii), opt.radii = 0.02; end

max_it=200;
[tar_vecs,tar_n,tar_lambda,tarData]=computeDescriper(tarData,opt);
[N,D]=size(tarData');
[M,D]=size(srcData');
if ~exist('sigma2','var') || isempty(sigma2) || (sigma2==0), 
    sigma2=(M*trace(tarData*tarData')+N*trace(srcData*srcData')-2*sum(srcData')*sum(tarData')')/(M*N*D);
end





iter=0;
WeightMatirx=[];

%%
while (iter<max_it) && (sigma2 > 1e-8) 
    TData=transform(srcData,opt.R,opt.t);
    [src_vecs,src_n,src_lambda,srcData]=computeDescriper(TData,opt);
    [M,D]=size(srcData');
    WeightMatrix=computeWeightMatrix(tar_vecs,src_vecs, tar_lambda,src_lambda,opt.beta);
    P_prior = WeightMatrix;
    iter=iter+1;
end
T=WeightMatrix;