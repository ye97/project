function [A]=reshape_NMD(A)
% [N,M,D]->[NM,D]
[N,M,D]=size(A);
A=permute(A,[2,1,3]);
A=reshape(A,[N*M,D]);
