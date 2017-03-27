function obj=Infill_Pseudo_EI(x, Kriging_model, f_min, point_added)
%--------------------------------------------------------------------------
% the pseudo EI criterion
%--------------------------------------------------------------------------
% Reference:
% D. Zhan, J. Qian, Y. Cheng, Pseudo expected improvement criterion for
% parallel, Journal of Global Optimization. doi:10.1007/s10898-016-0484-7
%--------------------------------------------------------------------------
% Zhan Dawei (zhandawei@hust.edu.cn)
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%--------------------------------------------------------------------------
% get the Kriging prediction and variance
if size(x,1)==1
    [yp,~,mse] = predictor(x,Kriging_model);
else
    [yp,mse] = predictor(x,Kriging_model);
end
s=sqrt(max(0,mse));
% calcuate the EI value
EI=(f_min-yp).*Gaussian_CDF((f_min-yp)./s)+s.*Gaussian_PDF((f_min-yp)./s);
%--------------------------------------------------------------------------
% if this is not the first infill point
if ~isempty(point_added) 
    % the scaling of x
    scaled_x = (x - Kriging_model.Ssc(1,:)) ./ Kriging_model.Ssc(2,:);
    scaled_point_added = (point_added - Kriging_model.Ssc(1,:)) ./ Kriging_model.Ssc(2,:);
    dx = scaled_x - scaled_point_added;
    correlation = feval(Kriging_model.corr, Kriging_model.theta, dx);
    % the Pseudo EI matrix
    EI=EI.*prod(1-correlation,1);
end
% the genetic algorithm tries to minimize the objective
obj=-EI;

end





