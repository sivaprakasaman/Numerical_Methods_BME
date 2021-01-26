% Andrew Sivaprakasam 
% BME 695- Numerical Methods
% Warm-up Assignment 1

%Iterate through each element i in a, swapping it with the element prior if element i-1 is greater than element i

function s = bubblesort(a)

n = length(a);
swapped = true;

while(swapped == true)
    swapped = false;
    for i = 2:(n)
        if(a(i-1)>a(i))
            first = a(i-1);
            second = a(i);
            
            a(i-1) = second;
            a(i) = first;
            
            swapped = true;
        end
        
    end
end 

s = a;

end

