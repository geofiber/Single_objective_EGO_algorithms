function [obj, con] = Fun_Pressure_Vessel(x)
% ----------------------------------------------------------------------------
% the pressure vessel design problem
% [1] Wilde, D. (1978) Globally Optimal Design,Wiley, NewYork
% [2] Rommel G. Regis. Constrained optimization by radial basis function
% interpolation for high-dimensional expensive black-box problems with
% infeasible initial points. Engineering Optimizaiton, 2014.
% fmin = 5804.45
% ----------------------------------------------------------------------------
x1 = x(:, 1); x2 = x(:, 2); x3 = x(:, 3); x4 = x(:, 4); 

obj = 0.6224*x1.*x3.*x4 + 1.7781*x2.*x3.^2 + 3.1661*x1.^2.*x4 + 19.84 * x1.^2.*x3;
con(:, 1) = -x1 + 0.0193*x3;
con(:, 2) = -x2 + 0.00954*x3;
con(:, 3) = plog(-pi*x3.^2.*x4 - 4*pi*x3.^3/3 + 1296000);

end