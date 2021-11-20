function [WeightMatrix] = computeWeightMatrix(tar_vecs,src_vecs, tar_lambda,src_lambda,beta)

WeightMatrix=[];
sizeX=size(tar_vecs);
sizeY=size(src_vecs);
sizeX=sizeX(1);
sizeY=sizeY(1);
pairs=0;
for i=1:sizeX
    for j=1:sizeY
        
          WeightMatrix{i,j}=norm(tar_vecs(i,:)-src_vecs(j,:))*norm(tar_lambda(i,:)-src_lambda(j,:));%局部向量做差距离*乘法向量夹角生成权重
%         if WeightMatrix{i,j}<0.15
%             pairs=pairs+1;
%         end
    end
end
%% 求和分配权重
WeightMatrix=cell2mat(WeightMatrix);
WeightMatrix=WeightMatrix./(-2*beta^2)
WeightMatrix=exp(WeightMatrix);
sumM=sum(WeightMatrix,2);
sumM==repmat(sumM,1,sizeY);
WeightMatrix= WeightMatrix./sumM;

