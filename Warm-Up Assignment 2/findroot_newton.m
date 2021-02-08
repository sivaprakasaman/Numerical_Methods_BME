function [root,x1,error] = findroot_newton(f,x0,epsilon)

syms x

f_prime = matlabFunction(diff(f(x)));
%f_prime2 = matlabFunction(diff(f_prime(x)));

i = 1;
x1(i) = x0;
error(i) = inf;
while(error>epsilon)
    i = i+1;
    
    x1(i) = x0 - f(x0)/f_prime(x0);
    error(i) = abs(f(x0)/(f_prime(x0)*x1(i)));
    x0 = x1(i);
end

root = x1;

end

