function [num_vari,design_space,optimum] = Test_Function(name)
%--------------------------------------------------------------------------
% single-objective optimization problem

switch name
    case 'Sphere'
        num_vari=6;
        design_space=[-5.12*ones(1,num_vari);5.12*ones(1,num_vari)];
        optimum=0;
        
    case 'SumSquare'
        num_vari=6;
        design_space=[-5.12*ones(1,num_vari);5.12*ones(1,num_vari)];
        optimum=0;
        
    case 'Rosenbrock'
        num_vari=6;
        design_space=[-5.12*ones(1,num_vari);5.12*ones(1,num_vari)];
        optimum=0;
        
    case 'Sixhump'
        num_vari=2;
        design_space=[-2,-2;2,2];
        optimum=-1.031628;
        
    case 'Branin'
        num_vari=2;
        design_space=[-5,0;10,15];
        optimum= 0.397887;
        
    case 'Sasena'
        num_vari=2;
        design_space=[0,0;5,5];
        optimum= -1.4565;
        
    case 'GoldPrice'
        num_vari=2;
        design_space=[-2,-2;2,2];
        optimum=3.0000;
        
    case 'Shekel5'
        num_vari=4;
        design_space=[0,0,0,0;10,10,10,10];
        optimum=-10.1532;
        
    case 'Shekel7'
        num_vari=4;
        design_space=[0,0,0,0;10,10,10,10];
        optimum=-10.4029;
        
    case 'Shekel10'
        num_vari=4;
        design_space=[0,0,0,0;10,10,10,10];
        optimum=-10.5364;
        
    case 'Hartman3'
        num_vari=3;
        design_space=[0,0,0;1,1,1];
        optimum=-3.8628;
        
    case 'Hartman6'
        num_vari=6;
        design_space=[0,0,0,0,0,0;1,1,1,1,1,1];
        optimum= -3.3224;
end
end




