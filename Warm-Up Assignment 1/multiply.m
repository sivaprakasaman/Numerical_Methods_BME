% Andrew Sivaprakasam 
% BME 695- Numerical Methods
% Warm-up Assignment 1

function product = multiply(x,y)
%Divide and Conquer aproach to int multiplication
%x and y in binary
%Testing

n = max(length(x),length(y))
x = dec2bin(bin2dec(x),n);
y = dec2bin(bin2dec(y),n);

    if n == 1
       product = bin2dec(x)*bin2dec(y);
    else
        m = n/2;
% 
%         x_R = x((ceil(m)+1):end)
%         x_L = x(1:floor(m))

        x_L = x(1:ceil(m))
        x_R = x((ceil(m)+1):n)
        

%         y_R = y((ceil(m)+1):end);
%         y_L = y(1:floor(m));
%         
                
        y_L = y(1:ceil(m))
        y_R = y((ceil(m)+1):n)
     
        P_1 = bin2dec(multiply(x_L, y_L))        
        %disp('P1_called')

        P_2 = bin2dec(multiply(x_R, y_R))
        %disp('P2_called')

        P_3 = bin2dec(multiply(dec2bin(bin2dec(x_L) + bin2dec(x_R),ceil(m)), dec2bin(bin2dec(y_L) + bin2dec(y_R),ceil(m))))
        %disp('P3_called')


        product = (P_1)*(2^(n)) + (P_3 - P_1 - P_2)*(2^(floor(m))) + P_2;

    end 
    product
    product = dec2bin(product);
    
end