% Andrew Sivaprakasam 
% BME 695- Numerical Methods
% Warm-up Assignment 1

%starting with element i, compare i against elements ahead. If there is a new minimum, swap it with element i

function s = selectionsort(a)

    k = length(a);
    
    for i=1:k-1
       tempMin = i;
       
       for j = (i+1):(k)
            if(a(j)<a(tempMin))
                tempMin = j;
            
            end 
       end
       
       if tempMin ~= i
            first = a(i);
            second = a(tempMin);
            
            a(i) = second;
            a(tempMin) = first;
            
       end 
    end
    s = a;
end
