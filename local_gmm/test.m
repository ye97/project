a=rand(3,3);
b=rand(2,3);
c=repmat(a,2,1).*repmat(b,3,1);
disp(a);
disp(b);
disp(c);
l=0;
for i=1:3
    for j=1:2
        l=l+1;
        d(l,:)=a(i,:).*b(j,:);
    end
end