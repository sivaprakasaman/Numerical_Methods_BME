%Andrew Sivaprakasam
%BME 695 | Numerical Methods
%Final Exam

clear all
close all

syms x y
%% Enter function here:

%need to save as two different function types fo rthis to work
% fx = @(x)(1.5-x(1)+x(1)*x(2))^2+(2.25-x(1)+x(1)*x(2)^2)^2 + (2.625 - x(1) + x(1)*x(2)^3)^2;;
% fxy = @(x,y) (1.5-x+x.*y)+(2.25-x+x.*y.^2).^2 + (2.625 - x + x.*y.^3).^2;


%mid is the start point (see last section of code)
%Rosenbrock N = 2
%  fx = @(x) 100.*(x(2)-x(1).^2).^2 + (1-x(1)).^2;
%  fxy = @(x,y) 100.*(y-x.^2).^2 + (1-x).^2;

%Sphere fxn
% fx = @(x) (x(1)-1)^2 + (x(2)+2)^2
% fxy = @(x,y) (x-1).^2 + (y+2).^2

% %Booth Fxn:
% fx = @(x) (x(1)+2*x(2)-7)^2 + (2*x(1) + x(2)-5)^2;
% fxy = @(x,y) (x+2.*y-7).^2 + (2.*x + y-5).^2;

% %Levi Fxn:
fx = @(x) sin(3*pi()*x(1))^2 + (1+sin(3*pi()*x(2))^2)*(x(1)-1)^2 + (1+sin(2*pi()*x(2))^2)*(x(2)-1)^2; 
fxy = @(x,y) sin(3.*pi().*x).^2 + (1+sin(3.*pi().*y).^2).*(x-1).^2 + (1+sin(2.*pi().*y).^2).*(y-1).^2; 

%% RSM Parameters:
bwidth = 0.01;
nbins = 10;
mid = [1,1];

%% Plotting the Surface

x = linspace(-6,6);
y = x;
[x,y] = meshgrid(x,y);

figure;
hold on
f_plot = feval(fxy,x,y);
surf(x,y,f_plot);
contour(x,y,f_plot);
xlabel('x')
ylabel('y')
zlabel('f(x,y)')
title('Function Surface/Contour Map n = 2')


%% Find True Min of fxn
opts = optimset('MaxFunEvals',20000,'MaxIter',20000);
[mins1,fval1] = fminsearch(fx,[4.5,4.5],opts);
[mins2,fval2] = fminsearch(fx,[-4.5,-4.5],opts);

%% Traverse and estimate min based on RSM

%Start with cc design centered at start point, find min
num_runs = 705;
%start point
%min_rsm = [4,4];
%true_val = [1,1];

%generate full factorial design
points = bwidth.*fullfact([nbins,nbins]);
points(:,1) = points(:,1)+(mid(1)-mean(points(:,1)));
points(:,2) = points(:,2)+(mid(2)-mean(points(:,2)));

func = feval(fxy,points(:,1),points(:,2));

m_error = 0;
minimum = max(func);
j = 1;

c = 1;

while(c>0.55)

points = bwidth.*fullfact([nbins,nbins]);
points(:,1) = points(:,1)+(mid(1)-mean(points(:,1)));
points(:,2) = points(:,2)+(mid(2)-mean(points(:,2)));    
    
func = feval(fxy,points(:,1),points(:,2));
coeffs = regress(func, horzcat(points(:,1), points(:,2), ones(length(points))));

estimates = coeffs(1)*points(:,1) + coeffs(2)*points(:,2) + coeffs(3);
actuals = feval(fxy,points(:,1),points(:,2));

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

    func = feval(fxy,points(:,1),points(:,2));
    
    polymodel = polyfitn(points, func,2);
    poly_s = polyn2sym(polymodel);
    
    poly_s2 = matlabFunction(poly_s);
    coeffs = polymodel.Coefficients;
    poly_s = @(x) coeffs(1)*x(1).^2 + coeffs(2)*x(1).*x(2) + coeffs(3)*x(1) + coeffs(4)*x(2).^2 + coeffs(5)*x(2) + coeffs(6);
    
    [min_rsm, fval_rsm(j)] = fmincon(poly_s,mid,[],[],[],[],[min(points(:,1)),min(points(:,2))],[max(points(:,1)),max(points(:,2))]);
    
    %in certain cases, the local minimum will be detected instead of the
    %global minimum...this means the search should be perturbed and re-run
    if((fval_rsm(j)>fval_rsm(j-1)))
        fval_rsm = fval_rsm(1:(j-1));
       %min_rsm_p(j,:) = min_rsm_p(:,1:j-1);
        break;
    end    
    
    min_rsm_p(j,:) = min_rsm;
    
    mid = min_rsm;
    
    %[min_rsm, fval_rsm(j)] = fmincon(poly_s,min_rsm,[],[],[],[],[-4.5 -4.5],[4.5 4.5]);    
    
    d = abs((fval_rsm(j)-fval_rsm(j-1))/fval_rsm(j));
    j = j+1;

end

    figure(1)
    c = linspace(1,10,num_runs);
    scatter3(min_rsm_p(:,1),min_rsm_p(:,2), fval_rsm(:),'filled');
    set(gca, 'CameraPosition', [-28.753761674980247,-58.83682331833516,941526.1673450879]);
    min_rsm