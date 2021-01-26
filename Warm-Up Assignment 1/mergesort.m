% Andrew Sivaprakasam 
% BME 695- Numerical Methods
% Warm-up Assignment 1

function s = mergesort(a)
    
    n = length(a);
    
    if n > 1        
        s = merge(mergesort(a(1:floor(n/2))),mergesort(a((floor(n/2)+1):n)));
    else
        s = a;
    end
    
end


function s = merge(arr1,arr2)

len1 = length(arr1);
len2 = length(arr2);

if isempty(arr1)
    s = arr2;
    return;
end

if isempty(arr2)
    s = arr1;
    return;
end

if arr1(1)<=arr2(1)
    m = merge(arr1(2:len1),arr2);
    s = [arr1(1), m];
else 
    m = merge(arr1,arr2(2:len2));
    s = [arr2(1), m];
end

end

