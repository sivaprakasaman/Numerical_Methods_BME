function product = multiply(xd,yd)
%Divide and Conquer aproach to int multiplication
%x and y in binary
%Testing

x = dec2bin(xd, 16);
y = dec2bin(yd, 16);

n = max(length(x),length(y));

    if n == 1
        product = str2num(x)*str2num(y);
    else
        m = n/2;

        x_L = bin2dec(x(ceil(m):end));
        x_R = bin2dec(x(1:floor(m)));

        y_L = bin2dec(y(ceil(m):end));
        y_R = bin2dec(y(1:floor(m)));

        P_1 = str2num(multiply(x_L, y_L))
        P_2 = str2num(multiply(x_R, y_R))
        P_3 = str2num(multiply(x_L + x_R, y_L + y_R))

        product = P_1*(2^n) + (P_3 - P_1 - P_2)*2^(m) + P_2; 
    end 
end

