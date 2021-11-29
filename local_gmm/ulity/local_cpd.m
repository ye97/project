function [T]=local_cpd(srcData,tarData,opt)
if ~isfield(opt,'H') || isempty(opt.H), opt.max_it = 10; end
if ~isfield(opt,'k') || isempty(opt.k), opt.k = 5; end
if ~isfield(opt,'downSample') || isempty(opt.downSample), opt.downSample=0; end
if ~isfield(opt,'s') || isempty(opt.s), opt.s = 1000; end
if ~isfield(opt,'radii') || isempty(opt.radii), opt.radii = 0.02; end
if ~isfield(opt,'radii') || isempty(opt.radii), opt.radii = 0.02; end
w=0.3;

[tar_n,tar_curvature,tar_localVec,tar_localDist]=findPointNormals(tarData,opt.k);
[N,D]=size(tarData);
[M,D]=size(srcData);
if ~exist('sigma2','var') || isempty(sigma2) || (sigma2==0), 
    sigma2=(M*trace(tarData*tarData')+N*trace(srcData*srcData')-2*sum(srcData)*sum(tarData)')/(M*N*D);
end

iter=0;
max_it=100;
WeightMatrix=[];
ntol=0;
tolerance=1e-9;
preE_con=0;

%% 
while (iter<max_it) %&& (ntol<1e-8) % && (sigma2 > 1e-8)% 
%     TData=srcData;
    TData=transform(srcData,opt.R,opt.t);
    [~,~,T_localVec,localDist]=findPointNormals(TData,opt.k);
    %% E step
%     paiMatrix=ones(N,M)./M;
    paiMatrix=computePai(tar_localVec,T_localVec, opt.beta);
    alpha=compute_alpha(tar_curvature,opt.alphamax);
    
    gloDist=compute_gloDist(tarData,TData,tar_n,alpha);
    P_prior =comput_Prior(paiMatrix,gloDist,w,sigma2);
    E_con=compute_E(P_prior,alpha,gloDist,sigma2);
    ntol=abs((E_con-preE_con)/E_con);
    preE_con=E_con;
    [opt.R,opt.t,sigma]=test_planeW(TData,tarData,tar_n,P_prior);
    sigma2=sigma^2;
   
    
    if(ntol <= tolerance)
        break;
    end
    %% M step
    iter=iter+1;
    disp(iter);
end
T=[opt.R,opt.t'];
cpd_plot_iter(tarData,TData);
hold on;
