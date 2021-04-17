%Andrew Sivaprakasam
%BME 695 | Numerical Methods
%Coding Assignment 2 Problem 1

function [x_out,y_out] = spline_2(x,y,x_out, end_cond)

if nargin == 3
    end_cond = 'natural';
end

n = length(x);
h = diff(x);

del_y = diff(y);

coeff_mat = zeros(n,n);
y_mat = zeros(n,1);

%set end conditions for various spline methods:

switch end_cond
    case 'natural'
        coeff_mat(1,:) = [1, zeros(1,n-1)];
        coeff_mat(end,:) = [zeros(1,n-1),1];
        y_mat(1) = 0;
        y_mat(end) = 0;
    case 'complete'
        coeff_mat(1,:) = [2*h(1),h(1), zeros(1,n-2)];
        coeff_mat(end,:) = [zeros(1,n-2),h(end-1),2*h(end-1)];
        
        %bad coding, only done for sake of assignment
        y1_prime = 1;
        yn_prime = 2;
        
        y_mat(1) = 3*del_y(1)/h(1) - 3*y1_prime;
        y_mat(end) = 3*yn_prime - 3*del_y(end)/h(end);
        
    case 'parabolic'
        coeff_mat(1,:) = [1, -1, zeros(1,n-2)];
        coeff_mat(end,:) = [zeros(1,n-2),1,-1];
        y_mat(1) = 0;
        y_mat(end) = 0;        
    
    case 'nak'
        coeff_mat(1,:) = [h(2), -(h(1)+h(2)), h(1), zeros(1,n-3)];
        coeff_mat(end,:) = [zeros(1,n-3),h(end-1),-(h(end-1)+ h(end-2)),h(n-1)];
        y_mat(1) = 0;
        y_mat(end) = 0;   
end

%Tridiagonal matrix
for i = 2:n-1
    coeff_mat(i,i-1:i+1) = [h(i-1),2*(h(i-1)+h(i)),h(i)];
    y_mat(i) = 3*(del_y(i)/h(i) - del_y(i-1)/h(i-1));
end

%Solve linear system of equations
c = inv(coeff_mat)*y_mat;

for j = 1:n-1
    d(j) = (c(j+1)-c(j))/(3*h(j));
    b(j) = del_y(j)/h(j) - c(j)*h(j) - d(j)*h(j)^2;
end

c = c(1:end-1)';

y_out = zeros(1,length(x_out));
for k = 1:n-1  
    y_out = y_out + (x_out>=x(k) & x_out<x(k+1)).*(y(k) + b(k).*(x_out-x(k)) + c(k).*(x_out-x(k)).^2 + d(k).*(x_out-x(k)).^3);    
end 
y_out(end) = y(end);
end

