function [alpha]=compute_alpha(tar_curvature,alpha_max)
lambda=1;
[N,D]=size(tar_curvature);
% alpha=((1-exp((3-1./tar_curvature).*lambda))./(1+exp((3-1./tar_curvature).*lambda))).*alpha_max;
alpha=ones(N,D);