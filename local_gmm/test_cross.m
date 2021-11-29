N=3;
M=4;
A=rand(N,3);
B=rand(M,3);
RESULT1=cross(repelem(A,M,1),repmat(B,N,1));
RESULT2=[];
for i=1:N
    for j=1:M
        RESULT2((i-1)*M+j,:)=cross(A(i,:),B(j,:));
    end
end
% % result=reshape(a,N*M,6);
% result=reshape(a,6,N*M);
% result=result';