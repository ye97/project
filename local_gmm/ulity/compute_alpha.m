function [alpha]=compute_alpha(tar_curvature,alpha_max)
lambda=1;
alpha=((1-exp((3-1./tar_curvature).*lambda))./(1+exp((3-1./tar_curvature).*lambda))).*alpha_max;
