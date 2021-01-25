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
