% Andrew Sivaprakasam
% Numerical Methods 
% Warm-Up Assignment 2

%Code generates all figures in the attached writeup

clear 
close all

%% Part 1.1| Bisection Method:

f1 = @(x)x.^3-0.165*x.^2+3.993e-4;
[root_f1_bsect, x1_root_b, error1_root_b] = findroot_bisect(f1,0.07,1,10e-6);

syms x
f2 = @(x)tan(x).*x.^(-2);
f2_prime = matlabFunction(diff(f2(x)));

[root_f2_bsect, x2_min_b, error2_min_b] = findroot_bisect(f2_prime,0.0001,2,10e-6);

%% Part 1.3 | Newton's Method:

[root_f1_newton,x1_root,error_root] = findroot_newton(f1,1,10e-6);
output_table = [(0:length(x1_root)-1)',x1_root,error_root];

[root_f2_newton,x1_minimization,error_minimization] = findroot_newton(f2_prime,1,10e-6);
output_table2 = [(0:length(x1_minimization)-1)',x1_minimization,error_minimization];
%% Comparing Convergence by Error:

hold on 
plot(error1_root_b,'LineWidth',1.5);
plot(error_root,'LineWidth',1.5);
hold off
title('Convergence of \epsilon in root finding')
ylabel('\epsilon');
xlabel('Iteration');
legend('Bisection','Newton');

figure
hold on 
plot(error2_min_b,'LineWidth',1.5);
plot(error_minimization,'LineWidth',1.5);
hold off
title('Convergence of \epsilon in minimization')
ylabel('\epsilon');
xlabel('Iteration');
legend('Bisection','Newton');


%% Part 1.5 | Task List 1 - C:

n_trials = 10;
epsilon = 10e-6:10e-2:1;
t_newton_f1 = zeros(n_trials, length(epsilon));
t_newton_f2 = zeros(n_trials, length(epsilon));
t_bisect_f1 = zeros(n_trials, length(epsilon));
t_bisect_f2 = zeros(n_trials, length(epsilon));

for i = 1:n_trials
   disp(['Iteration: ', num2str(i)]);
    for j = 1:length(epsilon)
        tic
            [~,~,~] = findroot_newton(f1,1,epsilon(j));
        t_newton_f1(i,j) = toc;
        
        tic
            [~,~,~] = findroot_newton(f2_prime,1,epsilon(j));
        t_newton_f2(i,j) = toc;
        
        tic
            [~,~,~] = findroot_bisect(f1,0.07,1,epsilon(j));
        t_bisect_f1(i,j) = toc;
        
        tic
            [~,~,~] = findroot_bisect(f2_prime,0.001,2,epsilon(j));
        t_bisect_f2(i,j) = toc;
    end
    
end

t_newton1_mean = mean(t_newton_f1);
t_newton1_std = std(t_newton_f1);
t_newton2_mean = mean(t_newton_f2);
t_newton2_std = std(t_newton_f2);

t_bisect1_mean = mean(t_bisect_f1);
t_bisect1_std = std(t_bisect_f1);

t_bisect2_mean = mean(t_bisect_f2);
t_bisect2_std = std(t_bisect_f2);

figure;
hold on 
errorbar(epsilon,t_newton1_mean,t_newton1_std,'LineWidth', 1.5);
errorbar(epsilon,t_bisect1_mean,t_bisect1_std,'LineWidth', 1.5);
hold off
title('Finding Zeros of x.^3-0.165*x.^2+3.993e^{-4}')
xlabel('\epsilon')
ylabel('Runtime (s)')
legend('Newton`s Method','Bisection Method');
ylim([-0.005,0.04]);
set(gca,'FontSize',11);

figure;
hold on 
errorbar(epsilon,t_newton2_mean,t_newton2_std,'LineWidth', 1.5);
errorbar(epsilon,t_bisect2_mean,t_bisect2_std,'LineWidth', 1.5);
hold off
title('Finding Minimum of tan(x).*x.^{-2}')
xlabel('\epsilon')
ylabel('Runtime (s)')
legend('Newton`s Method','Bisection Method');
ylim([-0.005,0.04]);
set(gca,'FontSize',11);

%% Part 2.1 | Euler's Method:

syms y t

y_prime1 = @(t,y)-y;
y_prime2 = @(t,y)2-exp(-4.*t)-2.*y;

[t_out_euler1,y_out_euler1] = euler_method(y_prime1,1,[0,5],0.1);
[t_out_euler2,y_out_euler2] = euler_method(y_prime2,1,[0,.5],0.1);

%% Part 2.1a | RK2 and RK4 

% RK2
[t_out_rk2_1,y_out_rk2_1] = rk2(y_prime1,1,[0,5],0.1,1,1); 
[t_out_rk2_2,y_out_rk2_2] = rk2(y_prime2,1,[0,.5],0.1,1,1);

% RK4
[t_out_rk4_1,y_out_rk4_1] = rk4(y_prime1,1,[0,5],0.1); 
[t_out_rk4_2,y_out_rk4_2] = rk4(y_prime2,1,[0,.5],0.1);

%% Analytical Solutions (for checking)
syms y(t)

fxn1 = diff(y,t) == -y;
cond = y(0) == 1;
y1Sol(t) = dsolve(fxn1,cond);

fxn2 = diff(y,t) == 2-exp(-4.*t)-2.*y;
cond = y(0) == 1;
y2Sol(t) = dsolve(fxn2,cond);

%% Plotting Comparisons:

figure;
hold on
plot(0:0.001:5,y1Sol(0:0.001:5),'k','LineWidth',1.5);
plot(t_out_euler1,y_out_euler1,'LineWidth',1.5);
plot(t_out_rk2_1,y_out_rk2_1,'LineWidth',1.5);
plot(t_out_rk4_1,y_out_rk4_1,'LineWidth',1.5);
hold off

title('Comparison of Solution to y` = -y');
xlabel('t');
ylabel('y');
legend('True','Euler','RK2','RK4');

figure;
hold on
plot(0:0.001:.5,y2Sol(0:0.001:.5),'k','LineWidth',1.5);
plot(t_out_euler2,y_out_euler2,'LineWidth',1.5);
plot(t_out_rk2_2,y_out_rk2_2,'LineWidth',1.5);
plot(t_out_rk4_2,y_out_rk4_2,'LineWidth',1.5);
hold off
title('Comparison of Solution to y` + 2y = 2-e^{-4t}');
xlabel('t');
ylabel('y');
legend('True','Euler','RK2','RK4');

%% Nelder-Mead:

%Minimize Rosenbrock's fxn using Nelder-Mead
f = @(x) 100.*(x(2)-x(1).^2).^2 + (1-x(1)).^2;
mins = fminsearch(f,[-1.2,1]);

