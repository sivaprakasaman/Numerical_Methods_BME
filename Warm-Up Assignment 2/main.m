% Andrew Sivaprakasam
% Numerical Methods 
% Warm-Up Assignment 2


%% Part 1.1| Bisection Method:

f1 = @(x)x.^3-0.165*x.^2+3.993e-4;
root_f1 = findroot_bisect(f1,0,0.07,0.001)

syms x
f2 = @(x)tan(x).*x.^(-2);
f2_prime = matlabFunction(diff(f2(x)));

root_f2_prime = findroot_bisect(f2_prime,0.00001,2,0.001);

%% Part 1.3 | Newton's Method:

[root_f1,x1_root,error_root] = findroot_newton(f1,0.0001,0.001);
output_table = [(0:length(x1_root)-1)',x1_root,error_root];

[root_f2,x1_minimization,error_minimization] = findroot_newton(f2_prime,0.0001,0.001);
output_table2 = [(0:length(x1_minimization)-1)',x1_minimization,error_minimization];
