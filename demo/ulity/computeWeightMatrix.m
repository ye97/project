function [WeightMatrix] = computeWeightMatrix(vecs_tar,vecs_src, tar_n,src_n,beta)

WeightMatrix=[];
sizeX=size(vecs_tar);
sizeY=size(vecs_src);
sizeX=sizeX(1);
sizeY=sizeY(1);
pairs=0;
for i=1:sizeX
    for j=1:sizeY
        
          WeightMatrix{i,j}=norm(vecs_tar{i,1}-vecs_src{j,1})*acos(dot(tar_n{i,1},src_n{j,1})/(norm(tar_n{i,1})*norm(src_n{j,1})));%局部向量做差距离*乘法向量夹角生成权重
          
          
          
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

