function s = insertionsort(a)

i = 2;

while(i<=length(a))
    j = i;
    
    while(j>1 && a(j-1)>a(j))
        
            first = a(j-1);
            second = a(j);
            
            a(j-1) = second;
            a(j) = first;
            
            j = j-1;
    end
    i = i+1;
end
    s = a;
end

