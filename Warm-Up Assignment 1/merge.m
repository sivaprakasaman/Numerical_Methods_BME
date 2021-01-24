function s = merge(arr1,arr2)

if isempty(arr1)
    s = arr2;
    return;
end

if isempty(arr2)
    s = arr1;
    return;
end

if arr1(1)<=arr2(1)
    m = merge(arr1(2:end),arr2);
    s = [arr1(1), m];
else 
    m = merge(arr1,arr2(2:end));
    s = [arr2(1), m];
end

end

