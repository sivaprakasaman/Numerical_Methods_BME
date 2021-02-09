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

[root_f2,x1_minimization,error_minimization] = findroot_newton(f2_prime,1,0.0000001);
output_table2 = [(0:length(x1_minimization)-1)',x1_minimization,error_minimization];

%% Part 2.1 | Euler's Method:

syms y t

y_prime1 = @(t,y)-y;
y_prime2 = @(t,y)2-exp(-4.*t)-2.*y;

[t_out_euler1,y_out_euler1] = euler_method(y_prime1,1,[0.1,5],0.1);
[t_out_euler2,y_out_euler2] = euler_method(y_prime2,1,[0.1,.5],0.1);

%% Part 2.1a | RK2 and RK4 

% RK2
[t_out_rk2_1,y_out_rk2_1] = rk2(y_prime1,1,[0.1,5],0.1,1,1); 
[t_out_rk2_2,y_out_rk2_2] = rk2(y_prime2,1,[0.1,.5],0.1,1,1);

% RK4
[t_out_rk4_1,y_out_rk4_1] = rk4(y_prime1,1,[0.1,5],0.1); 
[t_out_rk4_2,y_out_rk4_2] = rk4(y_prime2,1,[0.1,.5],0.1);

%% Analytical Solutions (for checking)
syms y(t)

fxn1 = diff(y,t) == -y;
cond = y(0) == 1;
y1Sol(t) = dsolve(fxn1,cond);

fxn2 = diff(y,t) == 2-exp(-4.*t)-2.*y;
cond = y(0) == 1;
y2Sol(t) = dsolve(fxn2,cond);

%% Plotting Comparisons:

%% Nelder-Mead:

%Minimize Rosenbrock's fxn using Nelder-Mead
f = @(x) 100.*(x(2)-x(1).^2).^2 + (1-x(1)).^2;
mins = fminsearch(f,[-1.2,1]);

