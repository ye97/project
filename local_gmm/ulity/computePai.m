function [paiMatrix] = computePai(tar_vecs,src_vecs,beta)


size_Tar=size(tar_vecs);
size_Src=size(src_vecs);
N=size_Tar(1);
M=size_Src(1);
D=size_Src(2);
% pairs=0;

tar_vecs=repmat(tar_vecs,1,M);
tar_vecs=reshape(tar_vecs,N,D,M);
tar_vecs=permute(tar_vecs,[1,3,2]);

src_vecs=repmat(src_vecs,1,N);
src_vecs=reshape(src_vecs,M,D,N);
src_vecs=permute(src_vecs,[3,1,2]);

% for i=1:sizeTar
%     for j=1:sizeSrc
%         
%           WeightMatrix{i,j}=norm(tar_vecs(i,:)-src_vecs(j,:))*norm(tar_lambda(i,:)-src_lambda(j,:));%局部向量做差距离*乘法向量夹角生成权重
% %         if WeightMatrix{i,j}<0.15
% %             pairs=pairs+1;
% %         end
%     end
% end

%% 计算相似度 ||tar_vecs-src_vecs||
paiMatrix=(tar_vecs-src_vecs).^2;
paiMatrix=sum(paiMatrix,3);
% patMatrix=sqrt(patMatrix);
%% 求和分配权重
% WeightMatrix=cell2mat(WeightMatrix);
paiMatrix=paiMatrix./(-2*beta^2);
paiMatrix=exp(paiMatrix);
sumM=sum(paiMatrix,2);
sumM=repmat(sumM,1,M);
paiMatrix= paiMatrix./sumM;

% 
% maxN=max(patMatrix);
% maxNM=max(maxN);

