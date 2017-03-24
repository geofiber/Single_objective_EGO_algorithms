%--------------------------------------------------------------------------
% The parallel Efficient Global Optimization (EGO) algorithm [1] using the
% pseudo Expected Improvement criterion (PEI) [2].
% the DACE toolbox of  Lophaven et al. (2002) [3] is used to build the kriging model
%--------------------------------------------------------------------------
% Reference:
% [1] Jones, D.R., Schonlau, M., Welch, W.J.: Efficient global optimization of
% expensive black-box functions. Journal of Global Optimization 13(4),
% 455-492 (1998).
%
% [2] D. Zhan, J. Qian, Y. Cheng, Pseudo expected improvement criterion for
% parallel, Journal of Global Optimization. doi:10.1007/s10898-016-0484-7
%
% [3]Lophaven SN, Nielsen HB, and Sodergaard J, DACE - A MATLAB Kriging
% Toolbox, Technical Report IMM-TR-2002-12, Informatics and Mathematical
% Modelling, Technical University of Denmark, 2002.
% Available at: http://www2.imm.dtu.dk/~hbn/dace/.
%--------------------------------------------------------------------------
% Zhan Dawei (zhandawei@hust.edu.cn)
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%--------------------------------------------------------------------------
clearvars;
addpath('dace');
addpath('test_problem')
%--------------------------------------------------------------------------
% you can choose from ¡®Sixhump¡¯, 'Branin', 'Sasena', 'GoldPrice',
% 'Shekel5', 'Shekel7', 'Shekel10', 'Hartman3', 'Hartman6', 'Sphere',
% 'SumSquare', 'Rosenbrock'
fun_name = 'GoldPrice';
% the number of initial design points
num_initial_sample = 50;
% the number of total allowed design points
num_total_sample = 100;
% the number of points selected in each iteration (cycle)
num_select = 4;
%--------------------------------------------------------------------------
% settings of the genetic algorithm
ga_population_size = 100;
ga_generation = 100;
ga_crossover_fraction = 0.9;
ga_option = gaoptimset( 'PopulationSize',ga_population_size, 'Generations',ga_generation,...
                                                   'StallGenLimit',ga_generation, 'CrossoverFraction',ga_crossover_fraction, 'Display', 'off');
    
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
Kriging_model = dacefit(sample_x,sample_y,'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
% the current iteration and evaluation
iteration = 0;
evaluation = size(sample_x,1);
% print the current information to the screen
fprintf(' iteration: %d, evaluation: %d, current best solution: %f: real optimum: %f\n', iteration, evaluation, f_min, optimum);

%--------------------------------------------------------------------------
% the iteration
while evaluation < num_total_sample
    %--------------------------------------------------------------------------
    % parameter settings for the optimization
    % how much candidates to be selected in this cycle
    num_candidate=min(num_select,num_total_sample-evaluation);
    % initial the candidate points and other parameters
    best_x = zeros(num_candidate,num_vari);
    best_EI = zeros(num_candidate, 1);
    point_added = [];
    
    %--------------------------------------------------------------------------
    % find the candidates based on pseudo EI criterion
    for ii = 1: num_candidate
        % the pseudo Expected Improvement criterion
        infill_criterion = @(x)pseudo_EI(x, Kriging_model, f_min, point_added);
        % find the point with the highest EI value using ga algorithm
        [best_x(ii,:),best_EI(ii,:)] = ga(infill_criterion,num_vari,[],[],[],[],design_space(1,:),design_space(2,:),[],ga_option);
        % update point_added
        point_added = best_x(1:ii,:);
    end
    % evaluating the candidate with the real function
    best_y = feval(fun_name,best_x);
    
    %--------------------------------------------------------------------------
    % add points and rebuild Kriging model
    % add the new point to design set
    sample_x = [sample_x;best_x];
    sample_y = [sample_y;best_y];
    % rebuild the Kriging model using the new design set
    Kriging_model = dacefit(sample_x,sample_y,'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
    % updating some parameters
    iteration = iteration+1;
    evaluation = size(sample_x,1);
    f_min = min(sample_y);
    % print the current information to the screen
    fprintf(' iteration: %d, evaluation: %d, current best solution: %f: real optimum: %f\n', iteration, evaluation, f_min, optimum);
    
end




