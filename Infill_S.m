function obj=Infill_S(x, Kriging_model)
%--------------------------------------------------------------------------
% the DACE toolbox of  Lophaven et al. (2002)  is used to predict value
%--------------------------------------------------------------------------
% REFERENCES
%Lophaven SN, Nielsen HB, and Sodergaard J, DACE - A MATLAB Kriging
%Toolbox, Technical Report IMM-TR-2002-12, Informatics and Mathematical
%Modelling, Technical University of Denmark, 2002.
%Available at: http://www2.imm.dtu.dk/~hbn/dace/.
%--------------------------------------------------------------------------
% get the Kriging prediction and variance
if size(x,1)==1
    [~,~,mse] = predictor(x,Kriging_model);
else
    [~,mse] = predictor(x,Kriging_model);
end
s=sqrt(max(0,mse));
% the genetic algorithm tries to minimize the objective
obj=-s;

end





