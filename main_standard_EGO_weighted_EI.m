%--------------------------------------------------------------------------
% the EGO algorithm using weighted EI criterion A. S¨®bester et al. (2005)
% is implemented
% the DACE toolbox of  Lophaven et al. (2002) [2] is used to build the kriging model
%--------------------------------------------------------------------------
% REFERENCES
 % [1] A. S¨®bester, S. Leary, A. Keane, On the Design of Optimization Strategies 
 % Based on Global Response Surface Approximation Models, Journal of Global 
 % Optimization 33(1) (2005) 31-59.
% [2] Lophaven SN, Nielsen HB, and Sodergaard J, DACE - A MATLAB Kriging
%Toolbox, Technical Report IMM-TR-2002-12, Informatics and Mathematical
%Modelling, Technical University of Denmark, 2002.
%Available at: http://www2.imm.dtu.dk/~hbn/dace/.
%--------------------------------------------------------------------------
% Zhan Dawei (zhandawei@hust.edu.cn)
% 2017.05.03, use the pso optimizer 
%--------------------------------------------------------------------------
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%--------------------------------------------------------------------------
clearvars;
addpath('dace');
addpath('test_problem')
%--------------------------------------------------------------------------
% setting of the problem
% you can choose from ¡®Sixhump¡¯, 'Branin', 'Sasena', 'GoldPrice',
% 'Shekel5', 'Shekel7', 'Shekel10', 'Hartman3', 'Hartman6', 'Sphere',
% 'SumSquare', 'Rosenbrock'
fun_name = 'GoldPrice';
% the number of initial design points
num_initial_sample = 20;
% the number of total allowed design points
num_total_sample = 100;
% the weighted in the weighted EI criterion 
weight_vector = [0.1, 0.3, 0.5, 0.7, 0.9];
%--------------------------------------------------------------------------
% this optimizer us the particle swarm optimization implemented in MATLAB 2016b
options=optimoptions('particleswarm','SwarmSize',100,'MaxIterations',100,'MaxStallIterations',100,'Display','off', 'UseVectorized', true );
% the genetic algorithm is repeated n times and the best solution is taken   
pso_repeat_time = 4;                                             
%--------------------------------------------------------------------------
 % get the information of the test problem                                           
[num_vari,design_space,optimum] = Test_Function(fun_name);
%--------------------------------------------------------------------------
% the 0th iteration
% initial design points using Latin hypercube sampling method
sample_x = repmat(design_space(1,:),num_initial_sample,1)+repmat(design_space(2,:)-design_space(1,:),num_initial_sample,1).*...
    lhsdesign(num_initial_sample,num_vari,'criterion','maximin','iterations',1000);
sample_y = feval(fun_name,sample_x);
% the current best solution
f_min = min(sample_y);
% build the initial Kriging model
kriging_model = dacefit(sample_x,sample_y,'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
% the current iteration and evaluation
iteration = 0;
evaluation = size(sample_x,1);
% print the current information to the screen
fprintf(' iteration: %d, evaluation: %d, current best solution: %f, real optimum: %f, weight : %f\n', iteration, evaluation, f_min, optimum,0);
%--------------------------------------------------------------------------
% the iteration
while evaluation < num_total_sample
    % the Expected Improvement criterion
    weight = weight_vector(iteration+1 - fix(iteration/5)*5);
    infill_criterion = @(x)Infill_Weighted_EI(x,kriging_model,f_min, weight);
    % find the point with the highest EI value using pso algorithm
    % and we run the pso optimizer pso_repeat_time times and take the best
    % result
    best_x_temp = zeros(pso_repeat_time, num_vari);
    best_EI_temp = zeros(pso_repeat_time , 1);
    for ii = 1 : pso_repeat_time
    [best_x_temp(ii,:), best_EI_temp(ii,:)] = particleswarm(infill_criterion, num_vari,design_space(1,:),design_space(2,:),options);       
    end
    [best_EI, ind] = max(-best_EI_temp);
    best_x = best_x_temp(ind,:);
    % evaluating the candidate with the real function
    best_y = feval(fun_name,best_x);    
    % add the new point to design set
    sample_x = [sample_x;best_x];
    sample_y = [sample_y;best_y];
    % rebuild the Kriging model using the new design set
    kriging_model = dacefit(sample_x,sample_y,'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
    % updating some parameters
    iteration = iteration+1;
    evaluation = size(sample_x,1);
    f_min = min(sample_y);
    % print the current information to the screen
    fprintf(' iteration: %d, evaluation: %d, current best solution: %f, real optimum: %f, weight: %f\n', iteration, evaluation, f_min, optimum,weight);
end




