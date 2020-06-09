function [obj,con] = Fun_Hesse(x)
% ----------------------------------------------------------------------------
% Hesse, R. 1973. A heuristic search procedure for estimating a global solution
% of nonconvex programming problems. Operations Research, 21:1267¨C1280.
% fmin = -310
% ----------------------------------------------------------------------------
x1 = x(:, 1); x2 = x(:, 2); x3 = x(:, 3); x4 = x(:, 4); x5 = x(:, 5); x6 = x(:, 6);

obj = -25*(x1 - 2).^2 - (x2 - 2).^2 - (x3 - 1).^2 - (x4 - 4).^2 - (x5 - 1).^2 - (x6 - 4).^2;
con(:,1) = (2 - x1 - x2)/2;
con(:,2) = (x1 + x2 - 6)/6;
con(:,3) = (-x1 + x2 - 2)/2;
con(:,4) = (x1 - 3*x2 -2)/2;
con(:,5) = (4 - (x3 - 3).^2 - x4)/4;
con(:,6) = (4 - (x5 - 3).^2 - x6)/4;


end