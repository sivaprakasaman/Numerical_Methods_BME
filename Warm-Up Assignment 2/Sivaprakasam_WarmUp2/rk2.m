% Andrew Sivaprakasam
% Numerical Methods 
% Warm-Up Assignment 2

function [t_out,y_out] = rk2(f,y0,tspan,h,a,b)

t_out = tspan(1):h:tspan(2);
y_out = zeros(length(t_out),1);
y_out(1) = y0;




for i = 2:length(t_out)
    
    k1 = feval(f,t_out(i-1),y_out(i-1));
    k2 = feval(f,t_out(i-1)+a*h,y_out(i-1)+b*h*k1);
    
    y_out(i) = y_out(i-1) + h*(.5*k1 + .5*k2);
    
end

end

