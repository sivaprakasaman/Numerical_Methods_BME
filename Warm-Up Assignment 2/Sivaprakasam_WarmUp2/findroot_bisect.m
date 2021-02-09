% Andrew Sivaprakasam
% Numerical Methods 
% Warm-Up Assignment 2

function [root,x_m2,error] = findroot_bisect(f,x_l,x_u,epsilon)

error = inf;
i = 0;
while(abs(error)>epsilon)
    i = i+1;
    x_m = (x_l+x_u)/2;
    prod = f(x_l)*f(x_m);

    if(prod==0)
        root = x_m;
        return;
    elseif(prod<0)
        x_u = x_m;    
    else
        x_l = x_m;
    end
    
    x_m2(i) = (x_l+x_u)/2;
    error(i) = abs((x_m2-x_m)/x_m2);
end

root = x_m2(i);
end

