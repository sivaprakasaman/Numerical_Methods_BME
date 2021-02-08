function root = findroot_bisect(f,x_l,x_u,epsilon)

error = inf;

while(abs(error)>epsilon)
    
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
    
    x_m2 = (x_l+x_u)/2;
    error = abs((x_m2-x_m)/x_m2);
end

root = x_m2;
end

