%Andrew Sivaprakasam
%Numerical Methods
%Warmup 1

function product = multiply2(x,y)

n = max(length(dec2bin(x)),length(dec2bin(y)));
n = 2^nextpow2(n);

x = dec2bin(x,n);
y = dec2bin(y,n);

if n == 1
    product = bin2dec(x)*bin2dec(y);
    return
else
    m = n/2;
    
    xL = bin2dec(x(1:ceil(m)));
    xR = bin2dec(x((ceil(m)+1):n));
    
    yL = bin2dec(y(1:ceil(m)));
    yR = bin2dec(y((ceil(m)+1):n));
    
    p1 = multiply2(xL, yL);
    p2 = multiply2(xR, yR);
    p3 = multiply2(xL+xR,yL+yR);
    
    product = p1*(2^(n)) + (p3 - p1 - p2)*(2^(floor(m))) + p2;
end     

end

