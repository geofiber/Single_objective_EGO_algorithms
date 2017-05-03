function obj=Infill_Weighted_EI(x, Kriging_model, f_min, weight)
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
[y,mse] = predictor(x,Kriging_model);
s=sqrt(max(0,mse));
% calcuate the EI value
EI1 = (f_min-y).*Gaussian_CDF((f_min-y)./s);
EI2 = s.*Gaussian_PDF((f_min-y)./s);
weighted_EI = weight * EI1 + (1-weight) * EI2;
% the genetic algorithm tries to minimize the objective
obj=- weighted_EI;

end





