% Andrew Sivaprakasam
% Numerical Methods 
% Warm-Up Assignment 2

function [t_out,y_out] = rk4(f,y0,tspan,h)

t_out = tspan(1):h:tspan(2);
y_out = zeros(length(t_out),1);
y_out(1) = y0;




for i = 2:length(t_out)
    
    k1 = h*feval(f,t_out(i-1),y_out(i-1));
    k2 = h*feval(f,t_out(i-1)+h/2,y_out(i-1)+k1/2);
    k3 = h*feval(f,t_out(i-1)+h/2,y_out(i-1)+k2/2);
    k4 = h*feval(f,t_out(i-1)+h,y_out(i-1)+k3);
    
    y_out(i) = y_out(i-1) + k1/6 + k2/3 + k3/3 + k4/6;
    
end

end
