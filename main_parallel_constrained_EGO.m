%--------------------------------------------------------------------------
% 1. the parallel constrained EGO algorithm [1] is implemented
%     for solving single-objective unconstrained expensive optimization problems
% 2. the DACE toolbox of  Lophaven et al. (2002) [2] is used to build the kriging model
% 3. The EI criterion is maximized by DE [3] algorithm.
%--------------------------------------------------------------------------
% REFERENCES
% [1] J. Qian, Y. Cheng, J. Zhang, J. Liu, and D. Zhan. A parallel constrained efficient
%     global optimization algorithm for expensive constrained optimization problems. 
%     Engineering Optimization, 2020, DOI:10.1080/0305215X.2020.1722118
% [2] Lophaven SN, Nielsen HB, and Sodergaard J, DACE - A MATLAB Kriging
%     Toolbox, Technical Report IMM-TR-2002-12, Informatics and Mathematical
%     Modelling, Technical University of Denmark, 2002.
%     Available at: http://www2.imm.dtu.dk/~hbn/dace/.
% [3] K. Price, R. M. Storn, and J. A. Lampinen, Differential evolution: 
%     a practical approach to global optimization: Springer Science & Business Media, 2006.
%     http://www.icsi.berkeley.edu/~storn/code.html
%--------------------------------------------------------------------------
% 2020.06.09
% zhandawei@swjtu{dot}edu{dot}cn
%--------------------------------------------------------------------------
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%--------------------------------------------------------------------------
clearvars;close all;
% settings of the problem
fun_name = 'Fun_Welded_Beam';
% get the information of the test problem
switch fun_name
    case 'Fun_Two_Member_Frame'
        num_vari = 3;    num_con = 2;  design_space = [2.5, 2.5, 0.1; 10, 10, 1];optimum = 703.916;
    case 'Fun_Gas_Transmission_Compressor'
        num_vari = 4;    num_con = 1;  design_space = [20,1,20,0.1;50,10,50,60]; optimum = 2964893.85;
    case 'Fun_Pressure_Vessel'
        num_vari = 4;    num_con = 3;  design_space = [0, 0, 0, 0; 1, 1, 50, 240];               optimum = 5804.45;
    case 'Fun_Welded_Beam'
        num_vari = 4;    num_con = 6;  design_space = [0.1, 0.1, 0.1, 0.1; 2.0, 10.0, 10.0, 2.0];   optimum = 1.725;
    case 'Fun_Hesse'
        num_vari = 6;     num_con = 6;  design_space = [0,0,1,0,1,0;5,4,5,6,5,10];                    optimum = -310;
    case 'Fun_Speed_Reducer'
        num_vari = 7;    num_con = 11;  design_space = [2.6, 0.7, 17, 7.3, 7.3, 2.9, 5.0; 3.6, 0.8, 28, 8.3, 8.3, 3.9, 5.5]; optimum = 2994.42;
    otherwise
        error('objective function is not defined!')
end
% the number of initial design points
num_initial = 20;
% the number of total allowed design points
max_evaluation = 100;
% the number of points selected in each iteration (cycle)
num_q = 5;
% the 0th iteration
% initial design points using Latin hypercube sampling method
sample_x = lhsdesign(num_initial, num_vari,'criterion','maximin','iteration',1000).*(design_space(2,:) - design_space(1,:)) + design_space(1,:) ;
[sample_y, sample_g] = feval(fun_name, sample_x);
% the number of total evaluations
evaluation = size(sample_x,1);
iteration = 0;
% record the f_min in each iteration
f_min = zeros(ceil((max_evaluation-num_initial)/num_q)+1 + 1,1);
% check is there is at least one feasible solution
index = sum(sample_g <= 0, 2) == num_con;
if sum(index) ~= 0
    f_min(1) = min(sample_y(index, :));
    fprintf('iteration: %d, evaluation: %d, best solution: %f, known optimum: %f\n', 0, evaluation, f_min(1), optimum);
else
    f_min(1) = 1E6;
    fprintf('iteration: %d, evaluation: %d, best solution: no feasiable solution, known optimum: %f\n', 0, evaluation, optimum);
end
%-------------------------------------------------------------------------
% beginning of the iteration
while evaluation < max_evaluation
    % if there is no feasiable solution, there is no need to build the
    % kriging model for the objective function (build anyway)
    kriging_con = cell(1, num_con);
    for ii = 1: num_con
        kriging_con{ii} = dacefit(sample_x, sample_g(:, ii),'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
    end
    kriging_obj = dacefit(sample_x,sample_y,'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
    % select the infill samples
    num_k = min(num_q,max_evaluation - evaluation);
    best_x = zeros(num_q, num_vari);
    point_added = [];
    for ii = 1 : num_q
        if sum(index) ~= 0
            infill_criterion = @(x)Infill_Pseudo_CEI(x, kriging_obj, kriging_con, f_min(iteration+1), point_added);
        else
            infill_criterion = @(x)Infill_Pseudo_PoF(x, kriging_con, point_added);
        end
        % find the candidate by maximizing the infill criterion
        [best_x(ii,:),best_EI] = DE(infill_criterion, num_vari, design_space(1,:), design_space(2,:), 100, 100);
        point_added = best_x(1:ii,:);
    end
    % evalaute the candidate points in parallel
    [best_y, best_g] = feval(fun_name,best_x);
    % add the new point to design set
    sample_x = [sample_x;best_x];
    sample_y = [sample_y;best_y];
    sample_g = [sample_g;best_g];
    % update the best solution
    % check is there is at least one feasible solution
    evaluation = size(sample_x,1);
    iteration = iteration + 1;
    index = sum(sample_g <= 0, 2) == num_con;
    if sum(index) ~= 0
        f_min(iteration+1) = min(sample_y(index, :));
        fprintf('iteration: %d, evaluation: %d, best solution: %f, known optimum: %f\n', iteration, evaluation, f_min(iteration+1), optimum);
    else
        f_min(iteration+1)= 1E6;
        fprintf(' iteration: %d, evaluation: %d, best solution: no feasiable solution, known optimum: %f\n',  iteration, evaluation, optimum);
    end
end





