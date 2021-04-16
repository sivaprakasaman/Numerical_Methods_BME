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
surf(x,y,f_beale)
figure;
contour(x,y,f_beale);

%% Rosenbrock Function

% syms x
% 
%Minimize Rosenbrock's fxn using Nelder-Mead
f = @(x,y) 100.*(y-x.^2).^2 + (1-x).^2;
mins = fminsearch(f,[-1.2,1]);

x = linspace(-100,100);
y = x;
[x,y] = meshgrid(x);

f_rosen = feval(f,x,y);
figure;
surf(x,y,f_rosen);
figure;
contour(x,y,f_rosen);

%% 2a Beale Minimization

syms x

opts = optimset('MaxFunEvals',20000,'MaxIter',20000);
f = @(x)(1.5-x(1)+x(1)*x(2))^2+(2.25-x(1)+x(1)*x(2)^2)^2 + (2.625 - x(1) + x(1)*x(2)^3)^2;
[mins1,fval1] = fminsearch(f,[4.5,4.5],opts);
[mins2,fval2] = fminsearch(f,[-4.5,-4.5],opts);

%% 2b Rosenbrock Minimization


