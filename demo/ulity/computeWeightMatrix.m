function [WeightMatrix] = computeWeightMatrix(Dx,Dy,beta)

WeightMatrix=[];
sizeX=size(Dx);
sizeY=size(Dy);
sizeX=sizeX(1);
sizeY=sizeY(1);
pairs=0;
for i=1:sizeX
    for j=1:sizeY
        
        WeightMatrix{i,j}=norm(Dx{i,1}-Dy{j,1})*norm(Dx{i,2}-Dy{j,2});
%         if WeightMatrix{i,j}<0.15
%             pairs=pairs+1;
%         end
    end
end
%% 按照指数分配权重
WeightMatrix=cell2mat(WeightMatrix);
WeightMatrix=WeightMatrix./(-2*beta^2)
WeightMatrix=exp(WeightMatrix);
sumM=sum(WeightMatrix,2);
sumM==repmat(sumM,1,sizeY);
WeightMatrix= WeightMatrix./sumM;

% %% 求和分配权重
% WeightMatrix=cell2mat(WeightMatrix);
% % WeightMatrix=WeightMatrix./(-2*beta^2)
% % WeightMatrix=exp(WeightMatrix);
% sumM=sum(WeightMatrix,2);
% sumM==repmat(sumM,1,sizeY);
% WeightMatrix= WeightMatrix./sumM;
