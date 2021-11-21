clc;
clear ;
close all;
load('p(Bun0).mat');
pc=15;
scannum=10;
res=10;
for i=1:scannum
    shape{1}=model{1,i}(1:end,:)
end


for i= 1:scannum
    mat=shape{i,1};
    ptCloud = pointCloud(mat);
    pcshow(ptCloud);
    pcwrite(ptCloud, strcat(num2str(i),".ply"));
end