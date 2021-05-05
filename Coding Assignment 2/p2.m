%Andrew Sivaprakasam
%BME 695 | Numerical Methods
%Coding Assignment 2 Problem 2

clear
close all;
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
figure(3);
hold on
surf(x,y,f_rosen);
contour(x,y,f_rosen);
set(gca, 'CameraPosition', [-50.293365694665845,-86.25634924965654,511952.0929456182]);
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
num_runs = 300;
%start point
%min_rsm = [4,4];
true_val = [1,1];

%Initial guess -> FF design
%bin-size

bwidth = 0.01;
nbins = 10;
mid = [4,4];

%generate full factorial design
points = bwidth.*fullfact([nbins,nbins]);
points(:,1) = points(:,1)+(mid(1)-mean(points(:,1)));
points(:,2) = points(:,2)+(mid(2)-mean(points(:,2)));


%f = @(x1,x2)(1.5-x1+x1.*x2).^2+(2.25-x1+x1.*x2.^2).^2 + (2.625 - x1 + x1.*x2.^3).^2;
f = @(x1,x2) 100.*(x2-x1.^2).^2 + (1-x1).^2;
func = feval(f,points(:,1),points(:,2));

m_error = 0;
minimum = max(func);
j = 1;

c = 1;

while(c>0.55)

points = bwidth.*fullfact([nbins,nbins]);
points(:,1) = points(:,1)+(mid(1)-mean(points(:,1)));
points(:,2) = points(:,2)+(mid(2)-mean(points(:,2)));

func = feval(f,points(:,1),points(:,2));
coeffs = regress(func, horzcat(points(:,1), points(:,2), ones(length(points))));

estimates = coeffs(1)*points(:,1) + coeffs(2)*points(:,2) + coeffs(3);
actuals = feval(f,points(:,1),points(:,2));

error = (estimates-actuals)./actuals;

%choose minimum based on least error of linear regression...if error is
%anything higher, the line needs to be shifted

c = corr(estimates,actuals);

[~,ind] = min(estimates);
minima = func(ind);

if(minima<minimum)
    mid = [points(ind,1), points(ind,2)];
    minimum = minima;
    min_rsm_p(j,:) = mid;
    fval_rsm(j) = minimum;
    j = j+1;
else
    break
end

end

clear points
d=1;
while(j <= num_runs && (d>0.001))

    points = ccdesign(2);
    points(:,1) = points(:,1)+(mid(1)-mean(points(:,1)));
    points(:,2) = points(:,2)+(mid(2)-mean(points(:,2)));

    func = feval(f,points(:,1),points(:,2));

    polymodel = polyfitn(points, func,2);
    poly_s = polyn2sym(polymodel);

    poly_s2 = matlabFunction(poly_s);
    coeffs = polymodel.Coefficients;
    poly_s = @(x) coeffs(1)*x(1).^2 + coeffs(2)*x(1).*x(2) + coeffs(3)*x(1) + coeffs(4)*x(2).^2 + coeffs(5)*x(2) + coeffs(6);

    [min_rsm, fval_rsm(j)] = fmincon(poly_s,mid,[],[],[],[],[min(points(:,1)),min(points(:,2))],[max(points(:,1)),max(points(:,2))]);

    if(fval_rsm(j)<fval_rsm(j-1))
        fval_rsm = fval_rsm(1:(j-1));
        break;
    end

    min_rsm_p(j,:) = min_rsm;

    %[min_rsm, fval_rsm(j)] = fmincon(poly_s,min_rsm,[],[],[],[],[-4.5 -4.5],[4.5 4.5]);

    d = abs(fval_rsm(j)-fval_rsm(j-1))/fval_rsm(j);


end

    func = feval(f,points(:,1),points(:,2));

    polymodel = polyfitn(points, func,2);
    poly_s = polyn2sym(polymodel);

    poly_s2 = matlabFunction(poly_s);
    coeffs = polymodel.Coefficients;
    poly_s = @(x) coeffs(1)*x(1).^2 + coeffs(2)*x(1).*x(2) + coeffs(3)*x(1) + coeffs(4)*x(2).^2 + coeffs(5)*x(2) + coeffs(6);

    [min_rsm, fval_rsm(j)] = fmincon(poly_s,min_rsm,[],[],[],[],[min(points(:,1)),min(points(:,2))],[max(points(:,1)),max(points(:,2))]);
    %[min_rsm, fval_rsm(j)] = fmincon(poly_s,min_rsm,[],[],[],[],[-4.5 -4.5],[4.5 4.5]);

    min_rsm_p(j,:) = min_rsm;

    j = j+1;
end

    figure(3)
    c = linspace(1,10,num_runs);
    scatter3(min_rsm_p(:,1),min_rsm_p(:,2), fval_rsm(:),'filled');
