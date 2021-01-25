% Andrew Sivaprakasam 
% BME 695- Numerical Methods
% Warm-up Assignment 1

function s = mergesort(a)
    
    n = length(a);
    
    if n > 1        
        s = merge(mergesort(a(1:floor(n/2))),mergesort(a((floor(n/2)+1):end)));
    else
        s = a;
    end
    
end


