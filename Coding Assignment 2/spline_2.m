function [x_out,y_out] = spline_2(x,y,x_out)

n = length(x);
h = diff(x);

del_y = diff(y);

%initial conditions for natural spline:
coeff_mat = zeros(n,n);
y_mat = zeros(n,1);
coeff_mat(1,:) = [1, zeros(1,n-1)];
coeff_mat(end,:) = [zeros(1,n-1),1];
y_mat(1) = 0;
y_mat(end) = 0;

for i = 2:n-1
    coeff_mat(i,i-1:i+1) = [h(i-1),2*(h(i-1)+h(i)),h(i)];
    y_mat(i) = 3*(del_y(i)/h(i) - del_y(i-1)/h(i-1));
end

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

