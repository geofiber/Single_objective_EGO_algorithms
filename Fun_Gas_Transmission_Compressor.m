function [obj,con] = Fun_Gas_Transmission_Compressor(x)
% ----------------------------------------------------------------------------
% Gas Transmission Compressor Design problem
% Beightler, C. S., D. T. Phillips, 1976. Applied Geometric Programming. Wiley, New York.
% global optimum is 2964893.85
% ----------------------------------------------------------------------------
x1 = x(:, 1); x2 = x(:, 2); x3 = x(:, 3); x4 = x(:, 4);

obj = 8.61E5*x1.^0.5.*x2.*x3.^(-2/3).*x4.^(-0.5) + 3.69E4*x3 + 7.72E8*x1.^(-1).*x2.^0.219 - 765.43E6*x1.^(-1);
con(:,1) = x4.*x2.^(-2) + x2.^(-2) -1;



end