%Andrew Sivaprakasam
%BME 695 | Numerical Methods
%Coding Assignment 3

%returns an error between observed/modeled.
function err = LVerr(obs_t, obs_N12, params, inits)

dt = 0.01;
tmin = 0; 
tmax = 20;

t = (tmin + dt):dt:tmax;

N1 = zeros(1,length(t));
N2 = zeros(1,length(t));

N1(1) = inits(1);
N2(1) = inits(2);

for i = 2:length(t)
    
    N1(i) = N1(i-1) + (params(1)*N1(i-1)-params(2)*N1(i-1)*(N2(i-1)))*dt;
    N2(i) = N2(i-1) + (params(3)*N1(i-1)*(N2(i-1))-params(4)*N2(i-1))*dt;

end

ind = obs_t/dt;

err = sum(([N1(ind);N2(ind)]-obs_N12).^2);

end

