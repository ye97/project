clc
clear
global Model Data
load bunny
%load p(Bun)
Model= bunny{1,1}'*1000;
Data= bunny{2,1}'*1000;
Model= Model(:,1:2:end);
Data= Data(:,1:2:end);
figure
plot3(Model(1,:),Model(2,:),Model(3,:),'.r');
hold on;
figure
plot3(Data(1,:),Data(2,:),Data(3,:),'.g');
R0 =eye(3);
t0= zeros(3,1);
tic
[R, t, ksi, mPhi] = FastTrICP(R0, t0, 100);
toc
Data2 = transform_to_global(Data, R, t);
figure 
plot3(Model(1,:),Model(2,:),Model(3,:),'.r');
hold on;
plot3(Data2(1,:),Data2(2,:),Data2(3,:),'.g');

