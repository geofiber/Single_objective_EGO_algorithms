function obj=Rosenbrock(x)
% Resenbrock test function


[m,n]=size(x);
obj=zeros(m,1);
for ii=1:m
    temp=0;
    for jj=1:n-1
        temp=temp+100*(x(ii,jj+1)-x(ii,jj)^2)^2+(x(ii,jj)-1)^2;
    end
    obj(ii,1)=temp;
end