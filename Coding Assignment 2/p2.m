%Andrew Sivaprakasam
%BME 695 | Numerical Methods
%Coding Assignment 2 Problem 2

%% Beale Function:


x = linspace(-4.5,4.5,100);
y = x;
[x,y] = meshgrid(x);
f_beale = (1.5-x+x.*y)+(2.25-x+x.*y.^2).^2 + (2.625 - x + x.*y.^3).^2;

% ff = f(xx,yy);
% levels = 100;
% contour(x,y,ff,levels);
figure;
hold on;
surf(x,y,f_beale)
contour(x,y,f_beale);
set(gca, 'CameraPosition', [-28.753761674980247,-58.83682331833516,941526.1673450879]);
hold off

xlabel('x')
ylabel('y')
zlabel('f(x,y)')
title('Beale Function Surface/Contour Map')
%% Rosenbrock Function

syms x
% 
%Minimize Rosenbrock's fxn using Nelder-Mead
f = @(x,y) 100.*(y-x.^2).^2 + (1-x).^2;
%mins = fminsearch(f,[-1.2,1]);

x = linspace(-6,6);
y = x;
[x,y] = meshgrid(x,y);

f_rosen = feval(f,x,y);
figure;
hold on
surf(x,y,f_rosen);
contour(x,y,f_rosen);
set(gca, 'CameraPosition', [-50.293365694665845,-86.25634924965654,511952.0929456182]);
hold off
xlabel('x')
ylabel('y')
zlabel('f(x,y)')
title('Rosenbrock Function Surface/Contour Map n = 2')
%% 2a Beale Minimization

syms x

opts = optimset('MaxFunEvals',20000,'MaxIter',20000);
f = @(x)(1.5-x(1)+x(1)*x(2))^2+(2.25-x(1)+x(1)*x(2)^2)^2 + (2.625 - x(1) + x(1)*x(2)^3)^2;

%Change fminsearch to fmincon specify constraints
[mins1,fval1] = fminsearch(f,[4.5,4.5],opts);
[mins2,fval2] = fminsearch(f,[-4.5,-4.5],opts);

%% 2b Rosenbrock Minimization for n = 2-4

f = @(x) (100.*(x(2)-x(1).^2).^2 + (1-x(1)).^2) + (100.*(x(3)-x(2).^2).^2 + (1-x(2)).^2) + (100.*(x(4)-x(3).^2).^2 + (1-x(3)).^2);
mins = fminsearch(f,[2,2,2,2]);

%% 3 Response Surface Optimization:

%Start with cc design centered at start point, find min

%Fit linear model

%If slope sign is negative, continue using linear model

%If slope sign changes to positive, switch to quadratic model

%Fit a quadratic model

%Minimize quadratic model using fminsearch

%stop after a given no of iterations
