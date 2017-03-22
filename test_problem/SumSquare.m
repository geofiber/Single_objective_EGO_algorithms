function obj=SumSquare(x)

A=repmat(1:1:size(x,2),size(x,1),1);
obj=sum(A.*x.^2,2);

end