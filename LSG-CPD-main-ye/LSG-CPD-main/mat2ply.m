clc;
clear ;
close all;
load('p(Dra0).mat');
pc=15;
for i= 1:pc 
    mat=model{1,i};
    ptCloud = pointCloud(mat);
    pcshow(ptCloud);
    pcwrite(ptCloud, strcat("dragon",num2str(i),".ply"));
end