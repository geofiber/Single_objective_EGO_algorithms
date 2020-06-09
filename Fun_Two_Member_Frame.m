function [obj, con] = Fun_Two_Member_Frame(x)
% ----------------------------------------------------------------------------
% design of a two-member frame
% Arora, J. S. (1989) Introduction to Optimum Design, McGraw-Hill Higher Education, NewYork.
% 2.5<=x1<=10
% 2.5<=x2<= 10
% 0.1<=x3<=1
% optimum found in above ref fmin = 703.916 at [7.798, 10, 0.1]
% ----------------------------------------------------------------------------
obj = zeros(size(x, 1) ,1);
con = zeros(size(x, 1) ,2);
for ii = 1 : size(x, 1)
    d = x(ii, 1); h = x(ii, 2); t = x(ii, 3);
    
    E = 3E7;
    G = 1.154E7;
    P = -10000;
    L = 100;
    
    A = (d - t)*(h - t);
    I = (d*h^3 - (d - 2*t)*(h - 2*t)^3)/12;
    J = 2*t*(d - t)^2*(h - t)^2/(d +h -2*t);
    
    K = [24, -6*L, 6*L;
        -6*L, 4*L^2 + G*J*L^2/(E*I), 0;
        6*L, 0, 4*L^2 + G*J*L^2/(E*I)];
    U = (K\[P;0;0])*L^3/(E*I);
    U1 =U(1,1);U2 =U(2,1);U3 =U(3,1);
    
    M1 = 2*E*I*(-3*U1 + U2*L)/L^2;
    M2 = 2*E*I*(-3*U1 + 2*U2*L)/L^2;
    T = -G*J*U3/L;
    
    sigma_1 = M1*h/(2*I);
    sigma_2 = M2*h/(2*I);
    tao = T/(2*A*t);
    
    obj(ii,1) = 2*L*(2*d*t + 2*h*t - 4*t^2);
    con(ii,1) = sqrt((sigma_1^2 + 3* tao^2))/40000 - 1;
    con(ii,2) = sqrt(sigma_2^2 + 3* tao^2)/40000 - 1;
    
end
