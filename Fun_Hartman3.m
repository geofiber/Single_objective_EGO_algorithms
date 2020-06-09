function obj = Fun_Hartman3(x)
%----------------------------------------------------------
% Hartman 3 Test Function for Nonlinear Optimization
% Taken from "Towards Global Optimisation 2",edited by L.C.W. Dixon and G.P.
% Szego, North-Holland Publishing Company, 1978. ISBN 0 444 85171 2
%
% 0 <= x1 <= 1
% 0 <= x2 <= 1
% 0 <= x3 <= 1
% fmin = -3.86278214782076
% xmin = [0.1, 0.55592003,0.85218259]
% http://www4.ncsu.edu/~definkel/research/index.html
%---------------------------------------------------------%

a=[3,10,30;
    0.1,10,35;
    3,10,30;
    0.1,10,35];
p=[0.3689,0.117,0.2673;
    0.4699,0.4387,0.747;
    0.1091,0.8732,0.5547;
    0.03815,0.5743,0.8828];
c=[1,1.2,3,3.2];
n=size(x,1);
obj=zeros(n,1);
for i=1:n
    obj(i,:)=-c*exp(-sum(a.*(repmat(x(i,:),4,1)-p).^2,2));
end
