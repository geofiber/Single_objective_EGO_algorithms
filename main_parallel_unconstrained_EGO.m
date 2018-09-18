%--------------------------------------------------------------------------
% The parallel Efficient Global Optimization (EGO) algorithm [1] using the
% pseudo Expected Improvement criterion (PEI) [2].
% for solving single-objective unconstrained expensive optimization problems
% the DACE toolbox of  Lophaven et al. (2002) [3] is used to build the kriging model
%--------------------------------------------------------------------------
% Reference:
% [1] Jones, D.R., Schonlau, M., Welch, W.J.: Efficient global optimization of
% expensive black-box functions. Journal of Global Optimization 13(4),
% 455-492 (1998).
% [2] D. Zhan, J. Qian, Y. Cheng, Pseudo expected improvement criterion for
% parallel. Journal of Global Optimization, 68(3):641-662.
% [3]Lophaven SN, Nielsen HB, and Sodergaard J, DACE - A MATLAB Kriging
% Toolbox, Technical Report IMM-TR-2002-12, Informatics and Mathematical
% Modelling, Technical University of Denmark, 2002.
% Available at: http://www2.imm.dtu.dk/~hbn/dace/.
%--------------------------------------------------------------------------
% Zhan Dawei (zhandawei@hust.edu.cn)
% 2017.05.03, use the pso optimizer
% 2018.09.18, run PSO optimizer only once for finding EI maximum
%--------------------------------------------------------------------------
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%--------------------------------------------------------------------------
clearvars;
%--------------------------------------------------------------------------
% setting of the problem
% you can choose from ¡®Sixhump¡¯, 'Branin', 'Sasena', 'GoldPrice',
% 'Shekel5', 'Shekel7', 'Shekel10', 'Hartman3', 'Hartman6'
fun_name = 'fun_Sixhump';
% get the information of the test problem
switch fun_name
    case 'fun_Sixhump'
        num_vari=2; design_space=[-2,-2;2,2];                  optimum=-1.031628;
    case 'fun_Branin'
        num_vari=2; design_space=[-5,0;10,15];                 optimum= 0.397887;
    case 'fun_Sasena'
        num_vari=2;  design_space=[0,0;5,5];                    optimum= -1.4565;
    case 'fun_GoldPrice'
        num_vari=2;  design_space=[-2,-2;2,2];                  optimum=3.0000;
    case 'fun_Shekel5'
        num_vari=4; design_space=[0,0,0,0;10,10,10,10];        optimum=-10.1532;
    case 'fun_Shekel7'
        num_vari=4; design_space=[0,0,0,0;10,10,10,10];         optimum=-10.4029;
    case 'fun_Shekel10'
        num_vari=4;   design_space=[0,0,0,0;10,10,10,10];       optimum=-10.5364;
    case 'fun_Hartman3'
        num_vari=3;   design_space=[0,0,0;1,1,1];               optimum=-3.8628;
    case 'fun_Hartman6'
        num_vari=6;  design_space=[0,0,0,0,0,0;1,1,1,1,1,1];    optimum= -3.3224;
    otherwise
        error('objective function is not defined!')
end
%--------------------------------------------------------------------------
% the number of initial design points
num_initial_sample = 20;
% the number of total allowed design points
num_total_sample = 100;
% the number of points selected in each iteration (cycle)
num_select = 5;
%--------------------------------------------------------------------------
% this optimizer us the particle swarm optimization implemented in MATLAB 2016b
options = optimoptions('particleswarm','SwarmSize',100,'MaxIterations',100,'MaxStallIterations',100,'Display','off', 'UseVectorized', true);
%--------------------------------------------------------------------------
% the 0th iteration
% initial design points using Latin hypercube sampling method
sample_x = design_space(1,:) + (design_space(2,:)-design_space(1,:)).*lhsdesign(num_initial_sample,num_vari,'criterion','maximin','iterations',1000);
sample_y = feval(fun_name,sample_x);
% the current best solution
f_min = min(sample_y);
% the current iteration and evaluation
iteration = 0;
evaluation = size(sample_x,1);
% print the current information to the screen
fprintf(' iteration: %d, evaluation: %d, current best solution: %f: real optimum: %f\n', iteration, evaluation, f_min, optimum);
%--------------------------------------------------------------------------
% the iteration
while evaluation < num_total_sample
    % build (or rebuild) the initial Kriging model
    kriging_model = dacefit(sample_x,sample_y,'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
    %--------------------------------------------------------------------------
    % parameter settings for the optimization
    % how much candidates to be selected in this cycle
    num_candidate = min(num_select,num_total_sample - evaluation);
    % initial the candidate points and other parameters
    best_x = zeros(num_candidate,num_vari);
    point_added = [];
    %--------------------------------------------------------------------------
    % find the candidates based on pseudo EI criterion
    for ii = 1: num_candidate
        % the pseudo Expected Improvement criterion
        infill_criterion = @(x)infill_pseudo_EI(x, kriging_model, f_min, point_added);
        % find the point with the highest EI value using ga algorithm
        best_x(ii,:)= particleswarm(infill_criterion, num_vari,design_space(1,:),design_space(2,:),options);
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
    % updating some parameters
    iteration = iteration+1;
    evaluation = size(sample_x,1);
    f_min = min(sample_y);
    % print the current information to the screen
    fprintf(' iteration: %d, evaluation: %d, current best solution: %f: real optimum: %f\n', iteration, evaluation, f_min, optimum);
end




