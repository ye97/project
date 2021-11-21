function [gloDist]=compute_Prior(tarData,srcData,,WeightMatirx, tar_n,tar_curvature)
N=size(tarData,1);
M=size(srcData,1);
D=size(srcData,2);

gloDist=N;