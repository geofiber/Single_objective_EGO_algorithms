function obj = Infill_Pseudo_CEI(x, kriging_obj, kriging_con, fmin, point_added)
%---------------------------------------------------
% the kriging prediction and varince
[u,mse] = predictor(x,kriging_obj);
s=sqrt(max(0,mse));
% the EI value
EI=(fmin-u).*Gaussian_CDF((fmin-u)./s)+s.*Gaussian_PDF((fmin-u)./s);
%---------------------------------------------------
% the number of constraints
num_con = length(kriging_con);
% the kriging prediction and varince
u_g = zeros(size(x,1), num_con);
mse_g = zeros(size(x,1), num_con);
for ii = 1: num_con
    [u_g(:, ii), mse_g(:, ii)] = predictor(x, kriging_con{ii});
end
s_g=sqrt(max(0,mse_g));
% the PoF value
PoF=prod(Gaussian_CDF((0-u_g)./s_g), 2);
%---------------------------------------------------
CEI = EI.*PoF;
%-----------------------------------
% if this is the first infill point
if ~isempty(point_added)
    % the scaling of x is the same for different objectives
    scaled_x = (x - kriging_obj.Ssc(1,:)) ./ kriging_obj.Ssc(2,:);
    scaled_point_added = (point_added - kriging_obj.Ssc(1,:)) ./ kriging_obj.Ssc(2,:);
    correlation = zeros(size(scaled_x,1),size(scaled_point_added,1));
    for ii =1:size(scaled_point_added,1)
        dx = scaled_x - scaled_point_added(ii,:);
        correlation(:,ii) = feval(kriging_obj.corr, kriging_obj.theta, dx);
    end
    % the Pseudo EI matrix
    CEI=CEI.*prod(1-correlation,2);
end
%--------------------------------------------------
% the objective is maximized
obj=-CEI;
end
