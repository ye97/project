
N=2;
M=3;
A=rand(N,3);
B=rand(M,3);
E1=cross(repelem(A,M,1),repmat(B,N,1));
E2=[];
for i=1:N
    for j=1:M
        E2(i,j,:)=cross(A(i,:),B(j,:));
    end
end
