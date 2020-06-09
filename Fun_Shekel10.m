function obj = Fun_Shekel10(x)
% this is the shekel function when m = 5.
% the design space is x_i = [0, 10], for i = 1,2,3,4
% the global optimum is x_opt = [4,4,4,4]
%                                             f_opt = -10.5364
a_matrix = [4, 4, 4, 4;
                        1, 1, 1, 1;
                        8, 8, 8, 8;
                        6, 6, 6, 6;
                        3, 7, 3, 7;
                        2, 9, 2, 9;
                        5, 5, 3, 3;
                        8, 1, 8, 1;
                        6, 2, 6, 2;
                        7, 3.6, 7, 3.6];
                    
c_vector = [0.1; 0.2; 0.2; 0.4; 0.4; 0.6; 0.3; 0.7; 0.5; 0.5];

obj = zeros(size(x,1),1);
for ii = 1: size(x,1)
obj(ii,:)=-sum((sum((repmat(x(ii,:),size(a_matrix,1),1)-a_matrix).^2,2)+c_vector).^-1);
end