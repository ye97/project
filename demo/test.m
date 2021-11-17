a=[1,2,3;
    4,5,6];
b= bsxfun(@rdivide,a,sum(a,2));
b=sum(a,2);


disp(b);
