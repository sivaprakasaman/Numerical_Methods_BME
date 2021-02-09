% Andrew Sivaprakasam
% Numerical Methods 
% Warm-Up Assignment 2

function [t_out,y_out] = euler_method(f,y0,tspan,h)

t_out = tspan(1):h:tspan(2);
y_out = zeros(length(t_out),1);
y_out(1) = y0;

    for i = 2:length(t_out)

        y_out(i) = y_out(i-1) + h*feval(f,t_out(i-1),y_out(i-1));

    end

end

