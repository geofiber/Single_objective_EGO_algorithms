%--------------------------------------------------------------------------
% 1. The parallel Efficient Global Optimization (EGO) algorithm [1] using the
%    pseudo Expected Improvement criterion (PEI) [2] for solving single-objective 
%    unconstrained expensive optimization problems
% 2. the DACE toolbox of  Lophaven et al. (2002) [3] is used to build the
%    kriging model.
% 3. The PEI criterion is maximized by DE [4] algorithm.
%--------------------------------------------------------------------------
% Reference:
% [1] Jones, D.R., Schonlau, M., Welch, W.J.: Efficient global optimization of
%     expensive black-box functions. Journal of Global Optimization 13(4),
%     455-492 (1998).
% [2] D. Zhan, J. Qian, Y. Cheng, Pseudo expected improvement criterion for
%     parallel. Journal of Global Optimization, 68(3):641-662.
% [3] Lophaven SN, Nielsen HB, and Sodergaard J, DACE - A MATLAB Kriging
%     Toolbox, Technical Report IMM-TR-2002-12, Informatics and Mathematical
%      Modelling, Technical University of Denmark, 2002.
%      Available at: http://www2.imm.dtu.dk/~hbn/dace/.
% [4] K. Price, R. M. Storn, and J. A. Lampinen, Differential evolution: 
%     a practical approach to global optimization: Springer Science & Business Media, 2006.
%     http://www.icsi.berkeley.edu/~storn/code.html
%--------------------------------------------------------------------------
% zhandawei@swjtu{dot}edu{dot}cn
% 2017.05.03, use the pso optimizer
% 2018.09.18, run PSO optimizer only once for finding EI maximum
% 2018.11.28  use DE optimizer for finding EI maximum
%--------------------------------------------------------------------------
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%--------------------------------------------------------------------------
clearvars; close all;
% setting of the problem
% you can choose from ¡®Sixhump¡¯, 'Branin', 'Sasena', 'GoldPrice',
% 'Shekel5', 'Shekel7', 'Shekel10', 'Hartman3', 'Hartman6'
fun_name = 'GoldPrice';
% get the information of the test problem
switch fun_name
    case 'Sixhump'
        num_vari=2; design_space=[-2,-2;2,2];                  optimum=-1.031628;
    case 'Branin'
        num_vari=2; design_space=[-5,0;10,15];                 optimum= 0.397887;
    case 'Sasena'
        num_vari=2;  design_space=[0,0;5,5];                    optimum= -1.4565;
    case 'GoldPrice'
        num_vari=2;  design_space=[-2,-2;2,2];                  optimum=3.0000;
    case 'Shekel5'
        num_vari=4; design_space=[0,0,0,0;10,10,10,10];        optimum=-10.1532;
    case 'Shekel7'
        num_vari=4; design_space=[0,0,0,0;10,10,10,10];         optimum=-10.4029;
    case 'Shekel10'
        num_vari=4;   design_space=[0,0,0,0;10,10,10,10];       optimum=-10.5364;
    case 'Hartman3'
        num_vari=3;   design_space=[0,0,0;1,1,1];               optimum=-3.8628;
    case 'Hartman6'
        num_vari=6;  design_space=[0,0,0,0,0,0;1,1,1,1,1,1];    optimum= -3.3224;
    otherwise
        error('objective function is not defined!')
end
%--------------------------------------------------------------------------
% the number of initial design points
num_initial = 20;
% the number of total allowed design points
max_evaluation = 200;
% the number of points selected in each iteration (cycle)
num_q = 5;
%--------------------------------------------------------------------------
% the 0th iteration
% initial design points using Latin hypercube sampling method
sample_x = repmat(design_space(1,:),num_initial,1) + repmat(design_space(2,:)-design_space(1,:),num_initial,1).*lhsdesign(num_initial,num_vari,'criterion','maximin','iterations',1000);
sample_y = feval(fun_name,sample_x);
% record the f_min in each iteration
f_min = zeros(ceil((max_evaluation-num_initial)/num_q)+1,1);
% the current best solution
f_min(1) = min(sample_y);
% the current iteration and evaluation
evaluation = size(sample_x,1);
iteration = 0;
% plot the iteration history
plot(iteration,f_min(1),'r-o');title(sprintf('iteration: %d, evaluations:%d',iteration,evaluation));drawnow;
% print the current information to the screen
fprintf(' iteration: %d, evaluation: %d, current best solution: %f: real optimum: %f\n', iteration, evaluation, f_min(1), optimum);
%--------------------------------------------------------------------------
% the iteration
while evaluation < max_evaluation 
    % build (or rebuild) the initial Kriging model
    kriging_model = dacefit(sample_x,sample_y,'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
    % initialize the candidate points and other parameters
    num_k = min(num_q,max_evaluation - evaluation);
    best_x = zeros(num_k,num_vari);
    point_added = [];
    % find the candidates based on pseudo EI criterion
    for ii = 1: num_k
        % the pseudo Expected Improvement criterion
        infill_criterion = @(x)Infill_Pseudo_EI(x, kriging_model, f_min(iteration+1), point_added);
        % find the point with the highest EI value using DE algorithm
        best_x(ii,:) = DE(infill_criterion, num_vari, design_space(1,:), design_space(2,:), 50, 200);
        % update point_added
        point_added = best_x(1:ii,:);
    end
    % evaluating the candidate with the real function
    best_y = feval(fun_name,best_x);
    % add the new point to design set
    sample_x = [sample_x;best_x];
    sample_y = [sample_y;best_y];
    % updating some parameters
    evaluation = size(sample_x,1);
    iteration = iteration + 1;
    f_min(iteration + 1) = min(sample_y);
    % plot the iteration history
    plot(0:iteration,f_min(1:iteration+1),'r-o');title(sprintf('iteration: %d, evaluations:%d',iteration,evaluation));drawnow;
    % print the current information to the screen
    fprintf(' iteration: %d, evaluation: %d, current best solution: %f: real optimum: %f\n', iteration, evaluation, f_min(iteration + 1), optimum);
end




