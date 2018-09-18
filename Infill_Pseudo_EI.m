function obj = infill_pseudo_EI(x, kriging_model, f_min, point_added)
%--------------------------------------------------------------------------
% the pseudo EI criterion
%--------------------------------------------------------------------------
% Reference:
% D. Zhan, J. Qian, Y. Cheng, Pseudo expected improvement criterion for
% parallel. Journal of Global Optimization, 68(3):641-662.
%--------------------------------------------------------------------------
% Zhan Dawei (zhandawei@hust.edu.cn)
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%--------------------------------------------------------------------------
% get the Kriging prediction and variance
[y,mse] = predictor(x,kriging_model);
s=sqrt(max(0,mse));
% calcuate the EI value
EI=(f_min-y).*gaussian_CDF((f_min-y)./s)+s.*gaussian_PDF((f_min-y)./s);
%--------------------------------------------------------------------------
% if this is not the first infill point
if ~isempty(point_added)
    % the scaling of x
    scaled_x = (x - kriging_model.Ssc(1,:)) ./ kriging_model.Ssc(2,:);
    scaled_point_added = (point_added - kriging_model.Ssc(1,:)) ./ kriging_model.Ssc(2,:);
    correlation = zeros(size(scaled_x,1),size(scaled_point_added,1));
    for ii =1:size(scaled_point_added,1)
        dx = scaled_x - scaled_point_added(ii,:);
        correlation(:,ii) = feval(kriging_model.corr, kriging_model.theta, dx);
    end
    % the Pseudo EI matrix
    EI = EI.*prod(1-correlation,2);
end
% the objective needs to be maximized
obj = -EI;

end





