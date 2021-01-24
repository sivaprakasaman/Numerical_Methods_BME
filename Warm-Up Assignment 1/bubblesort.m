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

