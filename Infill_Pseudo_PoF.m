function obj = Infill_Pseudo_PoF(x, kriging_con, point_added)
%------------------------------------
% the number of constraints
num_con = length(kriging_con);
% the kriging prediction and varince
u = zeros(size(x,1), num_con);
mse = zeros(size(x,1), num_con);
for ii = 1: num_con
    [u(:, ii), mse(:, ii)] = predictor(x, kriging_con{ii});
end
s=sqrt(max(0,mse));
%------------------------------------
% the PoF value
PoF=prod(Gaussian_CDF((0-u)./s), 2);
%-----------------------------------
% if this is the first infill point
if ~isempty(point_added)
    % the scaling of x
    scaled_x = (x - kriging_con{1}.Ssc(1,:)) ./ kriging_con{1}.Ssc(2,:);
    scaled_point_added = (point_added - kriging_con{1}.Ssc(1,:)) ./ kriging_con{1}.Ssc(2,:);
    correlation = zeros(size(scaled_x,1),size(scaled_point_added,1));
    for ii =1:size(scaled_point_added,1)
        dx = scaled_x - scaled_point_added(ii,:);
        correlation(:,ii) = feval(kriging_con{1}.corr, kriging_con{1}.theta, dx);
    end
    % the Pseudo EI matrix
    PoF=PoF.*prod(1-correlation,2);
end
%-----------------------------------
% the objective is maximized
obj=-PoF;
end
