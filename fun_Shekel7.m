function obj = fun_Shekel7(x)
% this is the shekel function when m = 5.
% the design space is x_i = [0, 10], for i = 1,2,3,4
% the global optimum is x_opt = [4,4,4,4]
%                                             f_opt = -10.4029
a_matrix = [4, 4, 4, 4;
                        1, 1, 1, 1;
                        8, 8, 8, 8;
                        6, 6, 6, 6;
                        3, 7, 3, 7;
                        2, 9, 2, 9;
                        5, 3, 5, 3];
                    
c_vector = [0.1; 0.2; 0.2; 0.4; 0.4; 0.6; 0.3];

obj = zeros(size(x,1),1);
for ii = 1: size(x,1)
obj(ii,:)=-sum((sum((x(ii,:)-a_matrix).^2,2)+c_vector).^-1);
end